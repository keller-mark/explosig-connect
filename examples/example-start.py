from explosig_connect import connect
import pandas as pd

conn = connect(session_id="4392af8f", empty=True, how=None, server_hostname='http://localhost:8100', client_hostname='http://localhost:8000')

df = pd.DataFrame(index=["SA 01", "SA 02", "SA 03"], columns=["Study", "Donor"], data=[
    {"Study": "Breast", "Donor": "DO 01"},
    {"Study": "Breast", "Donor": "DO 02"},
    {"Study": "Pancreatic", "Donor": "DO 02"}
])
conn.send_sample_metadata(df)

df = pd.DataFrame(index=["SA 01", "SA 02", "SA 03"], columns=["SBS", "DBS", "INDEL"], data=[
    {"SBS": 10, "DBS": 10, "INDEL": 10},
    {"SBS": 10, "DBS": 10, "INDEL": 1},
    {"SBS": 10, "DBS": 10, "INDEL": 1}
])
conn.send_mutation_type_counts(df)

df = pd.DataFrame(index=["COSMIC 01", "COSMIC 02", "COSMIC 03"], columns=["A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T"], data=[
    {"A[C>A]A": 0.2, "A[C>A]C": 0.3, "A[C>A]G": 0.1, "A[C>A]T": 0.4},
    {"A[C>A]A": 0.1, "A[C>A]C": 0.4, "A[C>A]G": 0.1, "A[C>A]T": 0.4},
    {"A[C>A]A": 0.25, "A[C>A]C": 0.25, "A[C>A]G": 0.25, "A[C>A]T": 0.25},
])
conn.send_signatures('SBS', df)

df = pd.DataFrame(index=["SA 01", "SA 02", "SA 03"], columns=["COSMIC 01", "COSMIC 02", "COSMIC 03"], data=[
    {"COSMIC 01": 200, "COSMIC 02": 300, "COSMIC 03": 100},
    {"COSMIC 01": 20, "COSMIC 02": 30, "COSMIC 03": 100},
    {"COSMIC 01": 100, "COSMIC 02": 300, "COSMIC 03": 200},
])
conn.send_exposures('SBS', df)