[![Build Status](https://travis-ci.org/keller-mark/explosig-connect.svg?branch=master)](https://travis-ci.org/keller-mark/explosig-connect)

## ExploSig Connect

Helpers for sending data from a Python environment to [ExploSig](https://github.com/lrgr/explosig) (via [ExploSig Server](https://github.com/lrgr/explosig-server)) for web-based interactive visualization.

### Installation

```sh
pip install explosig-connect
```

### Example Usage - Connecting to an existing session

```python
>>> from explosig_connect import connect

>>> # Connect using a session ID supplied by ExploSig.
>>> session_id = 'af6242f3'
>>> conn = connect(session_id)

>>> # Obtain the SBS mutation counts matrix associated with the session.
>>> sbs_counts_df = conn.get_mutation_type_counts('SBS')

>>> # You run some custom code to derive better signature exposures.
>>> exposures_df = my_exposures_computation_method(sbs_counts_df)

>>> # Send the new exposures back to ExploSig for visualization.
>>> conn.send_exposures('SBS', exposures_df)
```

### Example Usage - Starting a new session

```python
>>> from explosig_connect import connect

>>> # Start a new 'empty' session with no samples, signatures, etc. selected.
>>> # This will attempt to open ExploSig in a new browser tab that starts the session.
>>> conn = connect()

>>> # Send an SBS mutation counts matrix to visualize.
>>> conn.send_mutation_type_counts('SBS', my_sbs_counts_df)

>>> # Send an exposures matrix to visualize.
>>> conn.send_exposures('SBS', my_sbs_exposures_df)

>>> # Send a signatures matrix to visualize.
>>> conn.send_signatures('SBS', my_sbs_signatures_df)
```

### Development

Build and install from the current directory.

```sh
python setup.py sdist bdist_wheel
pip install .
```

Use the `hostname` parameter to specify a custom ExploSig instance.
There is also a `password` parameter to enable usage with protected instances.

```python
>>> conn = connect(session_id, hostname='http://localhost:8100')
```