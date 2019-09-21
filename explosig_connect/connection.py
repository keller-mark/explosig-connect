import requests
import json
import uuid
import webbrowser
import pandas as pd

class Connection:
    
    def post(self, data):
        payload = {
            'session_id': self.session_id,
            'token': self.token,
            'data': data,
        }
        r = requests.post(self.server_hostname + '/session-post', data=json.dumps(payload))
        r.raise_for_status()
        response = r.json()
        return response
    
    def send_sample_metadata(self, df):
        # df = samples x [ Study ]
        df.index = df.index.rename("sample_id")
        df = df.reset_index()
        df = df.rename(columns={"Study": "proj_id"})
        if "Donor" in df.columns.values.tolist():
            df = df.rename(columns={"Donor": "donor_id"})

        self.post({
            "data": {
                "data": {
                    "sample_meta": df.to_dict('records')
                },
                "scales": {
                    "sample_id": df["sample_id"].values.tolist(),
                    "proj_id": df["proj_id"].unique().tolist(),
                }
            }
        })
    
    def send_mutation_type_counts(self, df):
        # df = samples x [ SBS, DBS, INDEL ]
        df.index = df.index.rename("sample_id")
        mut_count_max = df.max().max()
        mut_count_sum_max = df.sum(axis=1).max()
        df = df.reset_index()

        self.post({
            "data": {
                "data": {
                    "mut_count": df.to_dict('records')
                },
                "scales": {
                    "mut_count": [0, int(mut_count_max)],
                    "mut_count_sum": [0, int(mut_count_sum_max)]
                }
            }
        })

    

class ConfigConnection(Connection):

    def __init__(self, session_id, token, server_hostname, client_hostname):
        self.session_id = session_id
        self.token = token
        self.server_hostname = server_hostname
        self.client_hostname = client_hostname

        self.config = self.get_config()

    def get_config(self):
        payload = {
            'session_id': self.session_id,
            'token': self.token,
        }
        r = requests.post(self.server_hostname + '/session-get', data=json.dumps(payload))
        r.raise_for_status()
        return json.loads(r.json()['state'])['config']
    
    def get_df(self, data_path, index_path, columns_path, index_col, mut_type=None, extra_config={}):
        payload = {
            'token': self.token,
            'projects': self.config['samples'],
            'tricounts_method': self.config['tricountsMethod'],
            'clinical_variables': self.config['clinicalVariables'],
            **extra_config,
        }
        # Try to provide signatures list even though using different naming conventions
        if mut_type != None:
            payload['mut_type'] = mut_type
            if mut_type == 'SBS':
                payload['signatures'] = self.config['signaturesSbs']
            elif mut_type == 'DBS':
                payload['signatures'] = self.config['signaturesDbs']
            elif mut_type == 'INDEL':
                payload['signatures'] = self.config['signaturesIndel']
        
        r_data = requests.post(self.server_hostname + data_path, data=json.dumps(payload))
        r_index = requests.post(self.server_hostname + index_path, data=json.dumps(payload))
        r_columns = requests.post(self.server_hostname + columns_path, data=json.dumps(payload))

        r_data.raise_for_status()
        r_index.raise_for_status()
        r_columns.raise_for_status()

        df = pd.DataFrame(data=r_data.json())
        df = df.set_index(index_col)

        index_df = pd.DataFrame(data=[], index=r_index.json(), columns=r_columns.json())
        df = df.reindex_like(index_df)

        return df
    
    def get_counts_by_category(self, mut_type):
        return self.get_df(
            '/plot-counts-by-category', 
            '/scale-samples', 
            '/scale-contexts', 
            'sample_id', 
            mut_type=mut_type,
        )

class EmptyConnection(Connection):

    def __init__(self, session_id, token, server_hostname, client_hostname):
        self.session_id = str(uuid.uuid4())[:8] if session_id == None else session_id[:8]
        self.token = token
        self.server_hostname = server_hostname
        self.client_hostname = client_hostname
    
    def open(self, how='auto'):
        if how == None:
            return
        assert(how in {'auto', 'nb_js', 'nb_link', 'browser'})

        url = self.client_hostname + '/#session-' + self.session_id

        if how in {'auto', 'browser'}:
            opened = webbrowser.open(url)
            if opened or how == 'browser':
                return
        try:
            from IPython import get_ipython
            if how in {'auto', 'nb_js'}:
                from IPython.display import display, Javascript
                js_block = "window.open('{}');".format(url)
                display(Javascript(js_block))
                return
            if how in {'auto', 'nb_link'}:
                from IPython.display import display, HTML
                html_block = "<a href='{}' target='_blank'>{}</a>".format(url, url)
                display(HTML(html_block))
                return
        except ImportError:
            print("Open the ExploSig session here: {}".format(url))
            return

