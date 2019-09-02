[![Build Status](https://travis-ci.org/keller-mark/explosig-connect.svg?branch=master)](https://travis-ci.org/keller-mark/explosig-connect)

## ExploSig Connect

Helpers for sending data from a Python environment to [ExploSig](https://github.com/lrgr/explosig) (via [ExploSig Server](https://github.com/lrgr/explosig-server)) for web-based interactive visualization.

### Installation

```sh
pip install explosig-connect
```

### Example Usage

```python
>>> from explosig_connect import connect

>>> # Connect using a session connection ID supplied by ExploSig.
>>> conn_id = 'af6242f3'
>>> conn = connect(conn_id)

>>> # Obtain the SBS mutation counts matrix associated with the session.
>>> sbs_counts_df = conn.get_counts_by_category('SBS')

>>> # You run some custom code to derive better signature exposures.
>>> exposures_df = my_exposures_computation_method(sbs_counts_df)

>>> # Send the new exposures back to ExploSig for visualization.
>>> conn.send_exposures('SBS', exposures_df)
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
>>> conn = connect(conn_id, hostname='http://localhost:8100')
```