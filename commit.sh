#!/bin/bash

# Automatically update timestamp before each commit

python version_update.py bard.py

git add .

git commit
