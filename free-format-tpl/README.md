# Free Format TPL files

The strict format of parameter files are prone to errors and conflicts. This project adds an option for users to develop parameters in free format then compile the parameters to current machine readable files.


## How does it work?
Free format parameters are put in folder maked as "free", such as param04-free and param08-free. Running "make" converts these files to corresponding folders param04 and param08. The old contents on param04 and param08 will be deleted.

Procedure:
1 prepare tpl file in param??-free folder (?? is the dielectric constant)
2 run "make param" to compile

File structure:
* `Makefile` - rules to compile toplogy files
* `compiletpl.py` - program to convert free format tpl to mcce readable tpl
* `param??-free/` - folders containing free format tpl files
* `param??/` - folders containing mcce readble tpl files, updated by "make" utility

## Implementation
Free format tpl files use json (JavaScript Object Notation) format.
