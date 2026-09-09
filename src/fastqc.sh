#!/bin/bash
gunzip *.gz
fastqc -o ../../results/ *.fastq