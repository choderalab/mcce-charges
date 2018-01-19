# Free Format TPL files

The strict format of parameter files are prone to errors and conflicts. This project adds an option for users
to develop parameters in free format, then compile the parameters to current machine readable files.


## How does it work?
Free format parameters are put in one folder and the converted parameter files will be put in another folder.

Procedure:
1 prepare a free format tpl file in a folder
2 run "tplconvert.py source_folder output_folder" to convert

## Implementation
Each parameter line consists of up to 3 key fields and one value. In free format, Keys and values are separated
by ":", and key fields are in quotation marks, separated by ",". In strict format, keys and value are solely
determined by position.

In addition, free format tpl files remove redundant entries like IATOM..

