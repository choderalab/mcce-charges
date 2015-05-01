# Free Format TPL files

The strict format of parameter files are prone to errors and conflicts. This project adds an option for users to develop parameters in free format, then compile the parameters to current machine readable files.


## How does it work?
Free format parameters are put in folder marked as "free", such as param04-free and param08-free. Running "make" converts these files to corresponding folders param04 and param08. The old contents on param04 and param08 will be deleted.

Procedure:
1 prepare tpl file in param??-free folder (?? is the dielectric constant)
2 run "make" under param??-free folder to compile.

File structure:
* `param??-free/` - folders containing free format tpl files
  * `Makefile` - rules to compile toplogy files
  * `compiletpl.py` - program to convert free format tpl to mcce readable tpl
* `param??/` - folders containing mcce readble tpl files, updated by "make" utility

## Implementation
Each parameter line consists of 3 keys and one value which may have sub fields. This is the same in free format tpl files.
Each key and value sub field is wrapped in quotation marks in free format tpl files, while they were solely determined
by position in strict format files.

In addition, free format tpl files remove redundant entries like IATOM.

