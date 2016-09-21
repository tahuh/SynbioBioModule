#!/usr/bin/python

"""
Custom memory usage check program
This code works for Linux Systen

Author : Sunghoon Heo


How to use

import os  ## This is required
from memory import MemoryUsageManager

pid = os.getpid()
manager = MemoryUsageManager(pid , unit = "kb") ## for Kilobytes. use Mb , Gb , Tb for other units
--- your runnign script here ---

manager.profile()

print manager.mem_usage ### This will give you memory in given units

mem_save_point1 = manager.save_point

--- your second code block here --

mamager.profile()

mem_save_point2 = manager.save_point

mem_diff = mem_save_point2 - mem_save_point1

print "Memory Difference : %f kb"%(mem_diff)

"""

import os
import re
import commands

class MemoryUsageManager(object) :
	def __init__(self , pid, unit="kb"):
		self.unit = unit
		self.denom = 0.0
		if unit == "kb" :
			self.denom = 1.0
		elif unit == "mb" or unit == "Mb" :
			self.denom = 1024.0
		elif unit == "gb" or unit == "Gb":
			self.denom = 1024.0 * 1024.0
		elif unit == "tb" or unit == "Tb":
			self.denom = 1024.0 * 1024.0 * 1024.0
		else:
			raise AttributeError("Invalid Memory unit. valid units are Kb, Mb, Gb, Tb")
		self.pid = pid

	def profile(self) :
		_proc_stat = "cat /proc/%s/status | grep VmSize"%(self.pid)
		o = commands.getoutput(_proc_stat)
		mem_usage = float(o.split()[1])
		self.mem_usage = str("%f %s"%(mem_usage / self.denom,self.unit))
		self.save_point = mem_usage / self.denom


