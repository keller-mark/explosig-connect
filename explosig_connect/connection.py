import requests
import json
import uuid
import webbrowser
import time
import pandas as pd
from pandas.api.types import is_numeric_dtype

class Connection:
    """
    Represents a connection to an ExploSig session,
    with functions for transforming and sending data.
    """
    
    def _post(self, data):
        # print(data) # for debugging 
        # print()
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
        """Send a dataframe containing sample metadata values to ExploSig.
        
        Parameters
        ----------
        df : `pandas.DataFrame`
            Dataframe with index of sample IDs. Columns are metadata variables.
            The following are recognized column names: {`Study`, `Donor`}.
        """

        df.index = df.index.rename("sample_id")
        df = df.reset_index()
        df = df.rename(columns={"Study": "proj_id"})
        if "Donor" in df.columns.values.tolist():
            df = df.rename(columns={"Donor": "donor_id"})

        self._post({
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
        """Send a dataframe containing mutation count values by mutation type to ExploSig.
        
        Parameters
        ----------
        df : `pandas.DataFrame`
            Dataframe with index of sample IDs. Columns are mutation types (`SBS`, `DBS`, `INDEL`).
        """
        # df = samples x [ SBS, DBS, INDEL ]
        df.index = df.index.rename("sample_id")
        mut_count_max = int(df.max().max())
        mut_count_sum_max = int(df.sum(axis=1).max())
        df = df.reset_index()

        self._post({
            "data": {
                "data": {
                    "mut_count": df.to_dict('records')
                },
                "scales": {
                    "mut_count": [0, mut_count_max],
                    "mut_count_sum": [0, mut_count_sum_max]
                }
            }
        })
    
    def send_signatures(self, mut_type, df, prob_max=None):
        """Send a dataframe containing signatures to ExploSig.
        
        Parameters
        ----------
        mut_type : `str`
            The mutation type corresponding to this set of signatures (`SBS`, `DBS`, `INDEL`).
        df : `pandas.DataFrame`
            Dataframe with index of signature names. Columns are mutation categories (`A[C>A]A`, etc.).
        prob_max : `None` or `'auto'`, optional
            How to compute the maximum y-value of signature plots. 
            If `None`, defaults to `0.2`. 
            If `'auto'`, set to the maximum value in the matrix.
            by default `None`
        """
        assert(mut_type in {'SBS', 'DBS', 'INDEL'})
        # df = signatures x categories
        if prob_max == 'auto':
            sig_prob_max = float(df.max().max())
        else:
            sig_prob_max = 0.2

        df.index = df.index.rename("sig_{}".format(mut_type))

        self._post({
            "data": {
                "scales": {
                    "cat_{}".format(mut_type): df.columns.values.tolist()
                }
            }
        })

        self._post({
            "data": {
                "scales": {
                    "sig_{}".format(mut_type): df.index.values.tolist()
                }
            }
        })

        self._post({
            "data": {
                "data": dict([
                    ("sig_{}_{}".format(mut_type, sig_name), df.loc[sig_name, :].to_frame(name="sig_prob_{}".format(mut_type)).reset_index().rename(columns={'index': "cat_{}".format(mut_type)}).to_dict('records')) for sig_name in df.index.values.tolist()
                ]),
                "scales": {
                    "sig_prob_{}".format(mut_type): [0.0, sig_prob_max]
                }
            }
        })
    
    def send_exposures(self, mut_type, df, send_sigs=False):
        """Send a dataframe containing exposures to ExploSig.
        
        Parameters
        ----------
        mut_type : `str`
            The mutation type corresponding to this set of signatures (`SBS`, `DBS`, `INDEL`).
        df : `pandas.DataFrame`
            Dataframe with index of sample IDs. Columns are signature names.
        send_sigs : `bool`, optional
            Whether to also send signature names with the exposures.
            Useful if not intending to call `send_signatures()`.
            by default `False`
        """
        assert(mut_type in {'SBS', 'DBS', 'INDEL'})

        if send_sigs:
            self._post({
                "data": {
                    "scales": {
                        "sig_{}".format(mut_type): df.columns.values.tolist(),
                    }
                }
            })
        
        mut_type = mut_type.lower() # Deal with inconsistent naming conventions
        df.index = df.index.rename("sample_id")
        exp_max = float(df.max().max())
        exp_sum_max = float(df.sum(axis=1).max())

        norm_df = pd.DataFrame(
            index=df.index.values.tolist(), 
            columns=df.columns.values.tolist(), 
            data=(df.values / df.values.sum(axis=1, keepdims=True))
        )
        norm_df.index = norm_df.index.rename("sample_id")
        exp_norm_max = float(norm_df.max().max())

        self._post({
            "data": {
                "data": {
                    "exposure_{}".format(mut_type): df.reset_index().to_dict('records'),
                    "exposure_{}_normalized".format(mut_type): norm_df.reset_index().to_dict('records'),
                },
                "scales": {
                    "exposure_{}".format(mut_type): [0.0, exp_max],
                    "exposure_sum_{}".format(mut_type): [0.0, exp_sum_max],
                    "exposure_{}_normalized".format(mut_type): [0.0, exp_norm_max],
                }
            }
        })
    
    def send_clinical_data(self, df, types={}, scales={}):
        """Send a dataframe containing clinical data.
        
        Parameters
        ----------
        df : `pandas.DataFrame`
            Dataframe with index of sample IDs. Columns are clinical variables.
        types : `dict`, optional
            A dict mapping column names to data types ('continuous' or 'categorical'), by default {}
            If a column name is not found in the dict, it is assumed that 
            numeric columns are continuous and string columns are categorical.
        scales : `dict`, optional
            A dict mapping column names to scale domains, by default {}
            If a column name is not found in the dict, it is assumed that 
            categorical column scales are simply a list of unique elements
            and continuous column scales are [min, max].
        """

        all_types = {}
        for col_name in df.columns.values.tolist():
            if col_name in types.keys():
                # check schema of types dict
                assert(types[col_name].lower() in {'continuous', 'categorical'})
                all_types[col_name] = types[col_name].lower()
            else:
                all_types[col_name] = 'continuous' if is_numeric_dtype(df[col_name]) else 'categorical'
        
        df.index = df.index.rename("sample_id")
  
        self._post({
            "data": {
                "data": {
                    "clinical_variable_type": [ {'variable': col_name, 'type': var_type} for col_name, var_type in all_types.items() ]
                },
                "scales": {
                    "clinical_variable": df.columns.values.tolist()
                }
            }
        })

        all_scales = {}
        
        for col_name, var_type in all_types.items():
            if col_name in scales.keys():
                # User provided the scale as a parameter
                all_scales[col_name] = scales[col_name]
            elif var_type == 'categorical':
                all_scales[col_name] = df[col_name].unique().tolist()
            elif var_type == 'continuous':
                all_scales[col_name] = [float(df[col_name].min()), float(df[col_name].max())]

        self._post({
            "data": {
                "data": {
                    "clinical_data": df.reset_index().to_dict('records')
                },
                "scales": all_scales
            }
        })
    
    def _send_gene_data(self, df, alteration_type):
        """Send a dataframe containing gene data for a specific alteration type.
        
        Parameters
        ----------
        df : `pandas.DataFrame`
            Dataframe with index of sample IDs. Columns are gene IDs.
        alteration_type : `str`
            One of {`mut_class`, `gene_expression`, `copy_number`}.
        """

        assert(alteration_type in {'mut_class', 'gene_expression', 'copy_number'})

        if alteration_type == 'mut_class':
            alteration_prefix = 'gene_mut'
        if alteration_type == 'gene_expression':
            alteration_prefix = 'gene_exp'
        if alteration_type == 'copy_number':
            alteration_prefix = 'gene_cna'

        df.index = df.index.rename("sample_id")

        genes = df.columns.values.tolist()

        self._post({
            "data": {
                "scales": {
                    alteration_prefix: genes,
                }
            }
        })

        for gene in genes:
            self._post({
                "data": {
                    "data": {
                        "{}_{}".format(alteration_prefix, gene): (df[gene].astype(str)
                            .to_frame().reset_index()
                            .rename(columns={ gene: alteration_type })
                            .to_dict('records'))
                    }
                }
            })
    
    def send_gene_mutation_data(self, df):
        """Send a dataframe containing gene mutation data.
        
        Parameters
        ----------
        df : `pandas.DataFrame`
            Dataframe with index of sample IDs. Columns are gene IDs.
        """
        self._send_gene_data(df, 'mut_class')
    
    def send_gene_expression_data(self, df):
        """Send a dataframe containing gene expression data.
        
        Parameters
        ----------
        df : `pandas.DataFrame`
            Dataframe with index of sample IDs. Columns are gene IDs.
        """
        self._send_gene_data(df, 'gene_expression')
    
    def send_copy_number_data(self, df):
        """Send a dataframe containing copy number data.
        
        Parameters
        ----------
        df : `pandas.DataFrame`
            Dataframe with index of sample IDs. Columns are gene IDs.
        """
        self._send_gene_data(df, 'copy_number')

    

class ConfigConnection(Connection):
    """
    Represents a connection to a previously-configured ExploSig session.
    """

    def __init__(self, session_id, token, server_hostname, client_hostname):
        self.session_id = session_id
        self.token = token
        self.server_hostname = server_hostname
        self.client_hostname = client_hostname

        self.config = self.get_config()

    def get_config(self):
        """Get the current data configuration as a `dict`.
        
        Returns
        -------
        `dict`
            A dictionary containing the selected samples, 
            signatures, clinical variables, and genes.
        """
        payload = {
            'session_id': self.session_id,
            'token': self.token,
        }
        r = requests.post(self.server_hostname + '/session-get', data=json.dumps(payload))
        r.raise_for_status()
        return json.loads(r.json()['state'])['config']
    
    def _get_df(self, data_path, index_path, columns_path, index_col, mut_type=None, extra_config={}):
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
    
    def get_mutation_type_counts(self):
        """Get the counts by mutation type dataframe associated with the current config.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and mutation types (SBS, DBS, INDEL) on the columns.
            Values are counts.
        """
        payload = {
            'token': self.token,
            'projects': self.config['samples']
        }

        data_path = '/plot-counts'
        r_data = requests.post(self.server_hostname + data_path, data=json.dumps(payload))
        r_data.raise_for_status()

        df = pd.DataFrame(data=r_data.json())
        df = df.set_index("sample_id")
        return df
    
    def get_mutation_category_counts(self, mut_type):
        """Get the counts by mutation category dataframe (for a particular mutation type) associated with the current config.
        
        Parameters
        ----------
        mut_type : `str`
            One of {`'SBS'`, `'DBS'`, `'INDEL'`}.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and mutation categories on the columns.
            Values are counts.
        """
        return self._get_df(
            '/plot-counts-by-category', 
            '/scale-samples', 
            '/scale-contexts', 
            'sample_id', 
            mut_type=mut_type,
        )
    
    def get_clinical_data(self):
        """Get the clinical data dataframe associated with the current config.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and clinical variables on the columns.
        """
        payload = {
            'token': self.token,
            'projects': self.config['samples']
        }
        data_path = '/plot-clinical'
        r_data = requests.post(self.server_hostname + data_path, data=json.dumps(payload))
        r_data.raise_for_status()

        df = pd.DataFrame(data=r_data.json())
        df = df.set_index("sample_id")

        return df
    
    def _get_gene_data(self, data_path, alteration_type):
        payload = {
            'token': self.token,
            'projects': self.config['samples']
        }
        index_path = '/scale-samples'
        genes = self.config['genes']

        r_index = requests.post(self.server_hostname + index_path, data=json.dumps(payload))
        r_index.raise_for_status()

        df = pd.DataFrame(data=[], index=r_index.json(), columns=genes)
        for gene_id in genes:
            r_data = requests.post(self.server_hostname + data_path, data=json.dumps({ **payload, 'gene_id': gene_id }))
            r_data.raise_for_status()

            gene_df = pd.DataFrame(data=r_data.json())
            gene_df = gene_df.set_index("sample_id")
            df[gene_id] = gene_df[alteration_type]

        return df
    
    def get_gene_mutation_data(self):
        """Get a dataframe containing mutation classes associated with the current config.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and genes on the columns.
            Values are mutation classes.
        """
        return self._get_gene_data('/plot-gene-mut-track', 'mut_class')
    
    def get_gene_expression_data(self):
        """Get a dataframe containing gene expression values associated with the current config.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and genes on the columns.
            Values are gene expression classes.
        """
        return self._get_gene_data('/plot-gene-exp-track', 'gene_expression')
    
    def get_copy_number_data(self):
        """Get a dataframe containing copy number values associated with the current config.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and genes on the columns.
            Values are copy number classes.
        """
        return self._get_gene_data('/plot-gene-cna-track', 'copy_number')
    
    def get_exposures(self, mut_type, tricounts_method=None):
        """Get the sample by signature exposures dataframe (for a particular mutation type) associated with the current config.
        
        Parameters
        ----------
        mut_type : `str`
            One of {`'SBS'`, `'DBS'`, `'INDEL'`}.
        tricounts_method: `str`, optional
            One of {`'By Study'`, `'None'`}. 
            Whether or not to normalize trinucleotides by frequency 
            (based on sequencing strategy of each selected cohort).
            By default, `'None'`.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and signature names on the columns.
            Values are counts (exposures).
        """
        assert(mut_type in {'SBS', 'DBS', 'INDEL'})

        payload = {
            'token': self.token,
            'projects': self.config['samples'],
            'mut_type': mut_type,
            'tricounts_method': (tricounts_method if tricounts_method != None else 'None'),
        }

        if mut_type == 'SBS':
            payload['signatures'] = self.config['signaturesSbs']
        elif mut_type == 'DBS':
            payload['signatures'] = self.config['signaturesDbs']
        elif mut_type == 'INDEL':
            payload['signatures'] = self.config['signaturesIndel']

        data_path = '/plot-exposures'
        r_data = requests.post(self.server_hostname + data_path, data=json.dumps(payload))
        r_data.raise_for_status()

        df = pd.DataFrame(data=r_data.json())
        df = df.set_index("sample_id")
        return df

class EmptyConnection(Connection):
    """
    Represents a connection to an "empty" ExploSig session.
    """

    def __init__(self, session_id, token, server_hostname, client_hostname):
        self.session_id = str(uuid.uuid4())[:8] if session_id == None else session_id[:8]
        self.token = token
        self.server_hostname = server_hostname
        self.client_hostname = client_hostname
    
    def open(self, how='auto'):
        """Attempts to open the session URL in a browser.
        Calls `webbrowser.open` if `how == 'browser'`.
        Outputs JavaScript if `how == 'nb_js'`.
        Outputs HTML if `how == 'nb_link'`.
        Otherwise, simply prints the URL.

        Parameters
        ----------
        how : `str`, optional
            One of {`'auto'`, `'nb_js'`, `'nb_link'`, `'browser'`}, by default `'auto'`
        """
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
    
    def _get_df(self, data_path, index_path, columns_path, index_col, extra_config={}):
        payload = {
            'token': self.token,
            **extra_config,
        }
        
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
    
    def get_mutation_type_counts(self, projects):
        """Get the counts by mutation type dataframe associated with the current config.

        Parameters
        ----------
        projects : `list` of `str`
            A list of sample cohort IDs.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and mutation types (SBS, DBS, INDEL) on the columns.
            Values are counts.
        """
        payload = {
            'token': self.token,
            'projects': projects
        }

        data_path = '/plot-counts'
        r_data = requests.post(self.server_hostname + data_path, data=json.dumps(payload))
        r_data.raise_for_status()

        df = pd.DataFrame(data=r_data.json())
        df = df.set_index("sample_id")
        return df
    
    def get_mutation_category_counts(self, mut_type, projects):
        """Get a mutation count dataframe (for a particular mutation type and set of sequencing projects).
        
        Parameters
        ----------
        mut_type : `str`
            One of {`'SBS'`, `'DBS'`, `'INDEL'`}.
        projects : `list` of `str`
            A list of sample cohort IDs.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and mutation categories on the columns.
            Values are counts.
        """
        extra_config = {
            'mut_type': mut_type,
            'projects': projects
        }
        return self._get_df(
            '/plot-counts-by-category', 
            '/scale-samples', 
            '/scale-contexts', 
            'sample_id', 
            extra_config=extra_config,
        )
    
    def get_clinical_data(self, projects):
        """Get a clinical data dataframe (for a particular set of sequencing projects).
        
        Parameters
        ----------
        projects : `list` of `str`
            A list of sample cohort IDs.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and clinical variables on the columns.
        """
        payload = {
            'token': self.token,
            'projects': projects
        }
        data_path = '/plot-clinical'
        r_data = requests.post(self.server_hostname + data_path, data=json.dumps(payload))
        r_data.raise_for_status()

        df = pd.DataFrame(data=r_data.json())
        df = df.set_index("sample_id")

        return df
    
    def _get_gene_data(self, genes, projects, data_path, alteration_type):
        payload = {
            'token': self.token,
            'projects': projects
        }
        index_path = '/scale-samples'

        r_index = requests.post(self.server_hostname + index_path, data=json.dumps(payload))
        r_index.raise_for_status()

        df = pd.DataFrame(data=[], index=r_index.json(), columns=genes)
        for gene_id in genes:
            r_data = requests.post(self.server_hostname + data_path, data=json.dumps({ **payload, 'gene_id': gene_id }))
            r_data.raise_for_status()

            gene_df = pd.DataFrame(data=r_data.json())
            gene_df = gene_df.set_index("sample_id")
            df[gene_id] = gene_df[alteration_type]

        return df
    
    def get_gene_mutation_data(self, genes, projects):
        """Get a dataframe containing mutation classes (for a particular set of genes and set of sequencing projects).
        
        Parameters
        ----------
        genes : `list` of `str`
            A list of gene IDs.
        projects : `list` of `str`
            A list of sample cohort IDs.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and genes on the columns.
            Values are mutation classes.
        """
        return self._get_gene_data(genes, projects, '/plot-gene-mut-track', 'mut_class')
    
    def get_gene_expression_data(self, genes, projects):
        """Get a dataframe containing gene expression values (for a particular set of genes and set of sequencing projects).
        
        Parameters
        ----------
        genes : `list` of `str`
            A list of gene IDs.
        projects : `list` of `str`
            A list of sample cohort IDs.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and genes on the columns.
            Values are gene expression classes.
        """
        return self._get_gene_data(genes, projects, '/plot-gene-exp-track', 'gene_expression')
    
    def get_copy_number_data(self, genes, projects):
        """Get a dataframe containing copy number values (for a particular set of genes and set of sequencing projects).
        
        Parameters
        ----------
        genes : `list` of `str`
            A list of gene IDs.
        projects : `list` of `str`
            A list of sample cohort IDs.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and genes on the columns.
            Values are copy number classes.
        """
        return self._get_gene_data(genes, projects, '/plot-gene-cna-track', 'copy_number')
    
    def get_exposures(self, projects, signatures, mut_type, tricounts_method=None):
        """Get the sample by signature exposures dataframe (for a particular mutation type) associated with the current config.
        
        Parameters
        ----------
        projects : `list` of `str`
            A list of sample cohort IDs.
        signatures : `list` of `str`
            A list of signature names.
        mut_type : `str`
            One of {`'SBS'`, `'DBS'`, `'INDEL'`}.
        tricounts_method: `str`, optional
            One of {`'By Study'`, `'None'`}. 
            Whether or not to normalize trinucleotides by frequency 
            (based on sequencing strategy of each selected cohort).
            By default, `'None'`.
        
        Returns
        -------
        `pandas.DataFrame`
            A dataframe with sample IDs on the index and signature names on the columns.
            Values are counts (exposures).
        """
        assert(mut_type in {'SBS', 'DBS', 'INDEL'})

        payload = {
            'token': self.token,
            'projects': projects,
            'mut_type': mut_type,
            'signatures': signatures,
            'tricounts_method': (tricounts_method if tricounts_method != None else 'None'),
        }

        data_path = '/plot-exposures'
        r_data = requests.post(self.server_hostname + data_path, data=json.dumps(payload))
        r_data.raise_for_status()

        df = pd.DataFrame(data=r_data.json())
        df = df.set_index("sample_id")
        return df

