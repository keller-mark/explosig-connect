import requests
import json
import pandas as pd

class Connection:

    def __init__(self, session_id, token, hostname):
        self.session_id = session_id
        self.token = token
        self.hostname = hostname

        self.config = self.get_config()

    def get_config(self):
        payload = {
            'session_id': self.session_id,
            'token': self.token,
        }
        r = requests.post(self.hostname + '/session-get', data=json.dumps(payload))
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
        
        r_data = requests.post(self.hostname + data_path, data=json.dumps(payload))
        r_index = requests.post(self.hostname + index_path, data=json.dumps(payload))
        r_columns = requests.post(self.hostname + columns_path, data=json.dumps(payload))

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

    def post(self, data):
        payload = {
            'session_id': self.session_id,
            'token': self.token,
            'data': data,
        }
        r = requests.post(self.hostname + '/session-post', data=json.dumps(payload))
        r.raise_for_status()
        response = r.json()
        return response

