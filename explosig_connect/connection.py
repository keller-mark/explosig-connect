import requests
import json
import pandas as pd

class Connection:

    def __init__(self, conn_id, token, hostname):
        self.conn_id = conn_id
        self.token = token
        self.hostname = hostname

        self.config = self.get_config()

    def get_config(self):
        payload = {
            'conn_id': self.conn_id,
            'token': self.token,
        }
        r = requests.post(self.hostname + '/session-get', data=json.dumps(payload))
        r.raise_for_status()
        return r.json()['config']
    
    def get_df(self, data_path, index_path, columns_path, index_col, extra_config={}):
        payload = {
            'token': self.token,
            **self.config,
            **extra_config,
        }
        r_data = requests.post(self.hostname + data_path, data=json.dumps(payload))
        r_index = requests.post(self.hostname + index_path, data=json.dumps(payload))
        r_columns = requests.post(self.hostname + columns_path, data=json.dumps(payload))

        r_data.raise_for_status()
        r_index.raise_for_status()
        r_columns.raise_for_status()

        full_df = pd.DataFrame(data=r_data.json())
        full_df = full_df.set_index(index_col)

        df = pd.DataFrame(data=[], index=r_index.json(), columns=r_index.columns())
        df = df.join(full_df, how='inner')

        return df
    
    def get_counts_by_category(self, mut_type):
        return self.get_df(
            '/plot-counts-by-category', 
            '/scale-samples', 
            '/scale-contexts', 
            'sample_id', 
            { 'mut_type': mut_type }
        )

    def post(self, data):
        payload = {
            'conn_id': self.conn_id,
            'token': self.token,
            'data': data,
        }
        r = requests.post(self.hostname + '/session-post', data=json.dumps(payload))
        r.raise_for_status()
        response = r.json()
        return response

