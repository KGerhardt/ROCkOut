import sys
import os
import subprocess

from .rocker_project_manager import project_manager
from .rocker_progress_tracker import progress_tracker

class run_pplacer:
	def __init__(self, dir = None):
		self.directory = dir
		self.reads = None
		self.ma = None
		
	def load_dir(self):
		rocker = project_manager(directory = project_directory, threads = threads)
		#This gets 
		rocker.parse_project_directory()
		#For cleaning, collect untagged reads
		rocker.parse_raw_reads()
		#get tagged reads, ready for alignment as QUERIES
		rocker.parse_tagged_reads()
		
	def generate_reference(self):
		pass
		