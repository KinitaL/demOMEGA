# demOMEGA

## Installation
### Cobra
1.Clone the repository
```shell
git clone
```
2. Create an environment
```shell
conda env create -f demOMEGA.yml
```
Note: if you use MacOS you have to specify the platform:
```shell
conda env create --platform osx-64 -f demOMEGA.yml
```
3. Activate the environment
```shell
conda activate demOMEGA
```

### Docker
#### From the repository
1. Clone the repository
```shell
git clone
```
2. Build and run the container
```shell
docker-compose run --rm --build demOMEGA
```
## Execution
Run synoSCO with default params (replace an example busco output with yours):
```shell
python3 demOMEGA.py
```

If we want to use nucleotide sequences it's enough to provide
only sequences or busco output. However, to use protein sequences it is a mandatory
to provide nucleotide sequences too (see config_example.yaml)


## Common errors
Problem:
```shell
bin/sh: shuf: command not found
```
Solution:
```shell
brew install coreutils
```



Written by: Nikita Leino

Aknowledgements. Many thanks to Giobbe Forni, whose work has pushed me to write this code.
