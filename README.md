# AMF analysis tutorial source

[![Copier](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/copier-org/copier/master/img/badge/badge-grayscale-inverted-border-orange.json)](https://github.com/copier-org/copier)

## Development instructions

The following sub-sections illustrate how to develop this documentation.

### Create the development environment

To build this documentation locally for development purposes, first create your development environment.

```
cd amf-tutorial
conda env create -n amf-tutorial --file environment-files/readthedocs.yml
conda activate amf-tutorial
q2doc refresh-cache
```



### Build the book

Next, build the book:

```
make html
```

(Alternatively, `make preview` or `make fast-preview` can speed up test builds.)

### Serve the book locally

Finally, run the following to serve the built documentation locally:

```
make serve
```

Have fun!
