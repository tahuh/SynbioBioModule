#!/usr/bin/python

"""
IlluminaSampleExtractor.py

This module only works for Synbio Lab index excel file format.

Author : Sunghoon Heo

Usage : python IlluminaSampleExtractor.py \
--excel1 <excel file that contains all index list> \
--excel2 <your sample list excel file> \
--out < output file name with path>

2016-12-02 Only MiniSeq avail
"""
try:
	import xlrd
except ImportError:
	raise ImportError("Python xlrd library required. Please tell administrator to intall.")

from argparse import ArgumentParser

class ExcelFiles(object):
	def __init__(self, excel1, excel2=None):
		self.generals = excel1
		if excel2 == None:
			self.specific = None
		self.specific = excel2


class ExcelReader(object):
	def __init__(self, EF):
		self.first = xlrd.open_workbook(EF.generals)
		if EF.specific == None : pass
		self.second = xlrd.open_workbook(EF.specific)
		self.content1 = []
		self.content2 = []
	def read_all(self):
		self.read_excel(self.first,1)
		if self.second != None: self.read_excel(self.second,2)
	
	def read_excel(self, obj,index=1):
		sheet = obj.sheet_by_index(0)
		numrows = sheet.nrows
		
		j = 0
		columns = [1,3,4]
		break_signal = False
		while j <  2:
			for r in xrange(0,numrows):
				### To solve index order we will first access by row and then access by columns
				### First columns required = 1,3,4 each represents p5, p7, sample name then 12, 14, 15
				data_pts = []
				for col in columns :
					try: cell_val = sheet.cell_value(r,col); data_pts.append(cell_val)
					except IndexError: break_signal=True; break
				if break_signal :
					break
				if index == 1 :
					self.content1.append(data_pts)
				elif index == 2 :
					self.content2.append(data_pts)
			if break_signal : break
			columns = [11,13,14]
			j += 1
			
class FormatterToRead(object):
	def __init__(self, reader,model='miniseq'):
		if model == None :
			raise AttributeError("Please specify model. miniseq, theragen avail currently(2016-12-02)")

#		self.reader = reader
#		self.reader2 = reader2
		self.excel_filled = self.simplify(reader,1)
		self.my_filled = self.simplify(reader,2)

	def miniseq_formatter(self, index) :
		index = index + 1
		prefix = str(index) + "_S" + str(index) +"_L001_R"
		fwd = prefix  + "1_001.fastq.gz"
		rev = prefix  + "2_001.fastq.gz"
		return fwd , rev
 	
	def simplify(self,R,index):
		if index == 1:
			extended = R.content1
		else :
			extended = R.content2
		new_list = []
		#print extended
		for entry in extended :
			
			e = entry[2]
			if e == None or e == '' or e == u''  or e == u'sample':
				continue
			else:
				new_list.append(entry)
		#print new_list
		return new_list

	def select(self, model='miniseq'):
		selected = []
		for data in self.my_filled:
			if data in self.excel_filled:
				i = self.excel_filled.index(data)
				if model == 'miniseq':
					fwd, rev = self.miniseq_formatter(i+1)
					selected.append((fwd, rev))
		# for i , data in enumerate(self.excel_filled):
			# if model == 'miniseq':
				# fwd, rev = self.miniseq_formatter(i)
				# if data in self.my_filled :
					# selected.append((fwd, rev))

		return selected


def main():
	parser = ArgumentParser()
	parser.add_argument("--excel1" , type=str,required=True)
	parser.add_argument("--excel2" , type=str,required=True)
	parser.add_argument("--out", type=str,required=True)
	args = parser.parse_args()
	excel_files = ExcelFiles(args.excel1, args.excel2)
	ExReader = ExcelReader(excel_files)
	ExReader.read_all()
	FormReader = FormatterToRead(ExReader)
	selected = FormReader.select()
	O = open(args.out,"w")
	for s in selected :
		
		O.write(s[0] + "\t" + s[1] + "\n")
	O.close()
	
if __name__ == "__main__":
	main()
