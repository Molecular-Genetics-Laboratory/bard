#!/bin/bash
# Convert XML layout definition into python code
# Need to attach all the signals and slots later
pyuic5 -x qtcreator_out.ui -o gui_template.py
