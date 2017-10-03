# fpr-analysis
Perform some analysis on the contents of the file provenance report. This probably won't work outside of OICR, but yay open science and open development.

Details are below for each script, but in short, to use them:

1. Generate the file provenance report (FPR) in JSON format. This can take a few minutes.
2. Run generate_csv_from_fpr.py on the FPR JSON file
3. Run boxplots_from_csv.py on the resulting CSV file

## File provenance report

Use the [file provenance report
client](https://github.com/oicr-gsi/pipedev/tree/master/pipedev-file-provenance-client)
to generate the JSON file. You will need a settings file set up to hit a
production environment. Grab this file from someone in pipeline.

    $ java -jar target/pipedev-file-provenance-client-2.2-jar-with-dependencies.jar --study DYS --study EACdysplasia --json --settings PRODUCTION.json -o fpr_dys.json

In this command, I'm only generating for DYS and EAC projects. These files can get pretty big. Running with --all turns into a 50G file.

## Generate CSV from FPR

I'm generating CSV primarily because the lab wants Excel spreadsheets to generate their own charts and graphs.

**Requirements**: python 2.7+, [ijson](https://pypi.python.org/pypi/ijson)

**Usage**: the script takes exactly one argument, the JSON FPR from the first step. The script must be run in a location with the appropriate file mounts to read the files listed in the FPR.

    python generate_csv_from_fpr.py FPR.json

**Output**: This will generate a file called FPR.json.csv in the current working directory with a header:

```
is_ffpe    aligned_bases   insert_stdev    soft_clip_bases reads_on_target average_read_length reads_per_start_point   insert_mean total_reads unmapped_reads  run_name    lane    barcode multiplex   target_size
```


## Generate Boxplots from CSV

Using the CSV from the last step, make boxplots comparing FFPE-Fresh Frozen for every column other than `is_ffpe` in the CSV, test for significance using an independent t-test.

**Requirements**: Python 2.7+, [NumPy](http://www.numpy.org/), [Matplotlib](https://matplotlib.org/), [SciPy](https://www.scipy.org/)


**Usage**: the script takes exactly one argument, the CSV from the second step. This step no longer needs FPR file access, but has many more dependencies.

    python boxplots_from_csv.py FPR.json.csv

**Output**: This will generate a PNG called FPR.json.csv.png in the current working directory with lots of box plots.


# License and help 

This project is distributed with the [MIT License](LICENSE) by Ontario Institute for Cancer Research and meee (Morgan Taschuk).

If you have any questions, please leave them on the bug tracker at https://github.com/morgantaschuk/fpr-analysis


