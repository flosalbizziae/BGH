#!/usr/bin/python3
#-*- coding: UTF-8 -*-  

import sys
import re
import itertools

#IUPAC Nucleotide ambiguity code
iupac={"A":"A",
"G":"G",
"C":"C",
"T":"T",
"Y":"[GA]",
"R":"[TC]",
"W":"[TA]",
"S":"[CG]",
"K":"[AC]",
"M":"[GT]",
"D":"[TCA]",
"V":"[TGC]",
"H":"[TGA]",
"B":"[GCA]",
"N":"[ATCG]",
"X":"[ATCG]",
"-":"-"}

#the function for linear digestion
def linear_digest(site,seq):
	#recognition site and sequence should be a string
	print(site,seq)
	fragments=[]
	if len(seq):
		for line in seq:
			previous=''
			#
			cut_start_pos=[m.start() for m in re.finditer(site,line)]
			cut_end_pos=[m.end() for m in re.finditer(site,line)]
			if len(cut_start_pos):
				print(cut_end_pos)#
				j=0
				for pos in cut_start_pos:
					digit_flg=0
					digit_len=0
					digit=0
					left_len=0
					right_len=0
					for x in cut_dict[site][0]:
						if x.isdigit():
							digit_flg+=1
							digit_len+=int(x)
					if digit_flg:
						left_len=len(cut_dict[site][0])-digit_flg*2+digit_len
						digit=1
						
					digit_flg=0
					digit_len=0
					for y in cut_dict[site][1]:
						if y.isdigit():
							digit_flg+=1
							digit_len+=int(y)				
					if digit_flg:
						right_len=len(cut_dcit[site][1])-digit_flg*2+digit_len
						digit=1
					
					if digit:
						#print(site)
						start_pos=pos+left_len
						end_pos=cut_end_pos[j]-right_len
					else:
						if (cut_end_pos[j]-pos+1)<(len(cut_dict[site][0])+len(cut_dict[site][1])):
							print('here')
							print(cut_end_pos[j] - pos)
							if j==0:
								start_pos=0
								end_pos=cut_start_pos[j+1]
							else:
								start_pos=cut_start_pos[j-1]
								end_pos=cut_start_pos[j+1]
							print(start_pos,end_pos)
						else:
							start_pos=pos+len(cut_dict[site][0])+1
							end_pos=start_pos
					#find the fragments according to the start and end position
					if j==0:
						if len(cut_start_pos)==1:
							fragments.append(line[:start_pos])
							fragments.append(line[end_pos:])
						else:
							fragments.append(line[:start_pos])
					else:
						if j==(len(cut_start_pos)-1):
							fragments.append(line[end_pos:start_pos])
							fragments.append(line[end_pos:])
						else:
							fragments.append(line[end_pos:start_pos])
					j+=1
			else:
				fragments+=[line]
				
		return fragments
#Check the number of input files
if len(sys.argv)<4:
	print("Your parameters are less then 2, please make sure to input both the restriction enzyme recognition site and the sequence file.")
	print("OPTIONS: -c	If this option is on, the program will digest the sequence as a circular DNA")
	print("OPTIONS: -l	If this option is on, the program will digest the sequence as a linear DNA")
	print("python digest.py [OPTIONS] <IN1 recognition site.txt> <IN2 sequence.fa>")
	print("INPUT1 EXAMPLE: \"BamHI\tC-GATCC|CCTAG-G\" One enzyme can recognize ambigurous sites and several enzyme is also allowed for one reaction.")
	print("INPUT2 EXMAPLE: refer to NCBI FASTA format specification")
else:
	#To avoid too much computation consumption caused by the combination of the enzymes, the enzyme number should be controlled of a reasonable range, here is 100
	fsite=open(sys.argv[2],"r")
	count=0
	while True:
		buffer = fsite.read(1000*1024)	#按照字节数读入文件内容给buffer（防止过大文件）
		if not buffer:
			break
		count+= buffer.count('\n')
	fsite.close( )
	if count>100:
		print("Your enzyme file contains more than 100 enzymes. Please make sure to input less than 100 enzymes.")
		sys.exit()
	else:
		#Deal the enzyme file and prepare the enzyme and sites dictionary
		fsite=open(sys.argv[2],"r")	#打开文件只读
		lines=fsite.readlines()	#注意内存的消耗
		lines=[line.strip().split("\t") for line in lines]	#列表解析实例-1
		enzymes=[x[0] for x in lines]	#列表解析实例-2
		sites=[x[1] for x in lines]	#列表解析实例-3
		tmp=zip(enzymes,sites)	
		enz_dict=dict(tmp)	#列表转字典
		fsite.close()	#使用后关闭文件
		
		site_list=[]
		endinfo=[]
		for site in enz_dict.values():
			tmp=(site.split("|"))[0]
			cut=tmp.split("-")	#t保留末端信息
			endinfo.append(cut)
			find="".join(str(i) for i in cut)
			tmp=find.replace("(","")
			find=tmp.replace(")","")
			#tiling short repeat sequences
			tmp=list(set(find))
			tmp.sort(key=find.index)
			core="".join(str(i) for i in tmp)
			if len(find)%len(core)==0 and len(core)<len(find):
				core="("+core+")+"	#考虑短重复序列的识别位点
				pattern=core
			else:
				core=find
				#treat the ambiguous sites
				pattern=''
				for i in core:
					if i in iupac:
						pattern+=iupac[i]	#处理iupac字符
					else:
						if i.isdigit():
							pattern+="{"+i+"}"	#处理多个字符匹配问题
						else:
							pattern+=i
			site_list.append(pattern)
		cut_dict=dict(zip(site_list,endinfo))
		
		#如果是多种酶需要排列组合，否则直接消化
		digest_order=[]
		#用递归实现所有元素排列的可能性
		digest_order=list(itertools.permutations(site_list))
		
		#print(cut_dict)
		val=''
		key=''
		if count < 100000:
			#小内存消耗（尤其是单条序列）
			fa_dict={}
			for line in open(sys.argv[3],"r"):	#逐行读入
				line=line.strip()	#去结尾换行
				regex = re.compile('^>')	#设定匹配模式
				m = re.match(regex,line)	#匹配
				if m is not None:
					key=line
					fa_dict[key]=''
				else:
					fa_dict[key]+=line
			#开始酶切
			result_dict={}
			for k in fa_dict:
				val=fa_dict[k].upper()
				if sys.argv[1] == "-l":
					for order in digest_order:
						Fragments=[val]
						for recognition_site in order:
							Fragments=linear_digest(recognition_site, Fragments)
						print("Fragments is :",Fragments)
						print("------------\n")
		else:
			#大内存消耗
			for line in open(sys.argv[3],"r"):	#逐行读入
				line=line.strip()	#去结尾换行
				regex = re.compile('^>')	#设定匹配模式
				m = re.match(regex,line)	#匹配并返回匹配结果
			
				if m is not None:
					if val is not None:
						#In this condition, the key and value are accordingly, and the sequence are merged to a complete string.
						#Dealing with the recognition site
						print(key)
						print(val)

					key=line
					val=''
				else:
					val+=line
		
		