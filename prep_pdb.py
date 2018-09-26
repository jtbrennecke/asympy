from pylab import *
from optparse import OptionParser

if __name__ == "__main__":
	"""The program renumbers the residues to start from 1 
	in any chain. The chains have to be labeled A,B,C,...
	"""
	parser = OptionParser()
	parser.add_option("-f","--file", dest="filename", help="pdb input file")
	parser.add_option("-o","--output", dest="output", help="pdb output file")
	(options,args) = parser.parse_args()

	chain = 'X'
	current_res = '0'
	resn = 0

	inf = open(options.filename,"r")
	outf = open(options.output,"w")
	for line in inf:
	    split = str.split(line)
	    if split[0] == "ATOM":
	
	        if chain != split[4]:
	            resn = 1
	            chain = split[4]
	            current_res = split[5]
	        else:
	            if current_res != split[5] :
	                resn += 1
	                current_res = split[5]
	
	        line_list = list(line)
	        format_str = "{:="+str(len(list(str(split[5]))))+"d}"
	        line_list[(26-len(list(str(split[5])))):26] = list(format_str.format(resn))
	        out_line=''.join(line_list)
	    else:
	        out_line = line
	    outf.write(out_line)
	inf.close()
	outf.close()
