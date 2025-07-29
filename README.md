[![logo](https://raw.githubusercontent.com/mistergroot/cladeaid/refs/heads/main/cladeaid_logo.png)](https://github.com/mistergroot/cladeaid)

Installation

```
    git clone https://github.com/mistergroot/cladeaid
    cd cladeaid
    conda env create -f environment.yml
    pip install -e .
```

Usage

```
    $cladeaid --help
    usage: cladeaid [-h] {lca,refinery} ...

    LCA assignment and downstream refinement tools

    positional arguments:
      {lca,refinery}
        lca           Run LCA assignment
        refinery      Run postprocessing/refinement

    options:
      -h, --help      show this help message and exit
```
