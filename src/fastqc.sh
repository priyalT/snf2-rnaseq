#!/bin/bash
gunzip *.gz
fastqc -o ../../results/ *.fq