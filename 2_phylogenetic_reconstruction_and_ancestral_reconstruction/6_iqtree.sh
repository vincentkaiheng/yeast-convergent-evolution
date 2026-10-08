#!/bin/bash

iqtree -s ./trimmed_sequences_sp_only -m MFP -T 30 -B 1000
