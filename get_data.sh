#7 mar 2023

wget https://raw.githubusercontent.com/cldf-datasets/phoible/master/cldf/values.csv
wget https://raw.githubusercontent.com/cldf-datasets/phoible/master/cldf/languages.csv

mkdir model_fits

python3 process_data.py
