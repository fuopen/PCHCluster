################################################ function part
################################################ define basic function

generate_evec<-function(prefix,bin_dir){
	#par.file<-paste0('par.',prefix)
	par.file<-gsub('(.*\\/)','\\1par.',prefix)
	cat(sprintf("genotypename:\t%s\n",paste0(prefix,'.geno')),file=par.file)
	cat(sprintf("snpname:\t%s\n",paste0(prefix,'.snp')),file=par.file,append=T)
	cat(sprintf("indivname:\t%s\n",paste0(prefix,'.ind')),file=par.file,append=T)
	cat(sprintf("evecoutname:\t%s\n",paste0(prefix,'.evec')),file=par.file,append=T)
	cat(sprintf("evaloutname:\t%s\n",paste0(prefix,'.eval')),file=par.file,append=T)
	cat("altnormstyle:\tNO\n",file=par.file,append=T)
	cat("numoutevec:\t10\n",file=par.file,append=T)
	cat("numoutlieriter:\t0",file=par.file,append=T)

	bin_dir<-gsub('\\/?$','/',bin_dir)
	smartpca<-paste0(bin_dir,'smartpca')
	twstats<-paste0(bin_dir,'twstats')
	cmd.smartpca<-paste0(smartpca,' -p ',par.file,' > ',prefix,'.log')
	cmd.twstats<-paste0(twstats,' -t ',bin_dir,'twtable -i ',prefix,'.eval -o ',prefix,'.tw')
	system(cmd.smartpca)
	system(cmd.twstats)
	file_evec<-paste(prefix,".evec",sep="")
	file_evec
}

read_evec<-function(file){
	dm<-read.table(file)
	data_matrix<-as.matrix(dm[,c(2:3)])
	data_matrix
}

get_p<-function(file){
	a<-read.table(file)
	pvalue=a[1,5]
	pvalue
}

correct1_firstPC<-function(PC1){
	PC1_ind1<-which(PC1<=0)
	PC1_ind2<-which(PC1>0)
	PC1_order=order(PC1)
	PC1_order_neg=PC1[PC1_order[1:length(PC1_ind1)]]
	PC1_order_pos=PC1[PC1_order[(length(PC1_ind1)+1):length(PC1)]]
	if(length(PC1_ind1)>length(PC1_ind2)){
#print("situation1")
		var_vector=c()
		for(i in 1:(length(PC1_ind2)-3)){
			var_vector=c(var_vector,(var(PC1_order_neg)+var(PC1_order_pos)))
			PC1_order_neg=c(PC1_order_neg,PC1_order_pos[1])
	        PC1_order_pos=PC1_order_pos[-1]
		}
		num_extra_ind1=which.min(var_vector)-1
			if(num_extra_ind1>0){
				num_ind1=length(PC1_ind1)+num_extra_ind1
				PC1_ind1=PC1_order[1:num_ind1]
				PC1_ind2=PC1_order[-c(1:num_ind1)]
			}
	}
	if(length(PC1_ind1)<length(PC1_ind2)){
#print("situation2")
		var_vector=c()
		for(i in 1:(length(PC1_ind1)-3)){
			var_vector=c(var_vector,(var(PC1_order_neg)+var(PC1_order_pos)))
			PC1_order_neg=PC1_order_neg[-length(PC1_order_neg)]
			PC1_order_pos=c(PC1_order_neg[length(PC1_order_neg)],PC1_order_pos)
		}
		num_extra_ind2=which.min(var_vector)-1
		print(num_extra_ind2)
		if(num_extra_ind2>0){
			num_ind2=length(PC1_ind2)+num_extra_ind2
			print(num_ind2)
			PC1_ind1=PC1_order[1:(length(PC1_order)-num_ind2)]
			PC1_ind2=PC1_order[-c(1:(length(PC1_order)-num_ind2))]
		}
	}
	return(list(PC1_ind1=PC1_ind1,PC1_ind2=PC1_ind2))	
}


sep_file<-function(prefix,data_matrix,if_ini=FALSE){
	PC1<-data_matrix[,1]
	inds<-correct1_firstPC(PC1)
	pc1_1_ind<-sort(inds[[1]])
	pc1_2_ind<-sort(inds[[2]])
	samplesize<-length(PC1)	

	generate_ind(pc1_1_ind,prefix,1,if_ini)
	generate_ind(pc1_2_ind,prefix,2,if_ini)

	generate_gen(pc1_1_ind,prefix,samplesize,1,if_ini)
	generate_gen(pc1_2_ind,prefix,samplesize,2,if_ini)

	generate_snp(prefix,1,if_ini)
	generate_snp(prefix,2,if_ini)
}

generate_ind<-function(pc_ind,prefix,i,if_ini=FALSE){
	ind<-read.table(paste(prefix,".ind",sep=""))
	ind_sep<-ind[pc_ind,]
	ind_sep[,2]<-"F"	
	if(if_ini==TRUE){
		prefix=paste(prefix,"_",sep="")
	}
	if(i==1){
		write.table(ind_sep,file=paste(prefix,"0.ind",sep=""),quote=F,row.names=F,col.names=F)
		pops<-ind_sep[,3]
		write.table(unique(pops),file=paste(prefix,"0.pop",sep=""),quote=F,col.names=F)
	}
	else{
		write.table(ind_sep,file=paste(prefix,"1.ind",sep=""),quote=F,row.names=F,col.names=F)
		pops<-ind_sep[,3]
		write.table(unique(pops),file=paste(prefix,"1.pop",sep=""),quote=F,col.names=F)
	}
}

generate_gen<-function(pc_ind,prefix,samplesize,i,if_ini=FALSE){
	gen<-read.fwf(paste(prefix,".geno",sep=""),width=rep(1,samplesize))
	gen_sep<-gen[,pc_ind]
	if(if_ini==TRUE){
		prefix=paste(prefix,"_",sep="")
	}
	if(i==1){
		write.table(gen_sep,file=paste(prefix,"0.geno",sep=""),quote=F,row.names=F,col.names=F,sep="")
	}
	else{
		write.table(gen_sep,file=paste(prefix,"1.geno",sep=""),quote=F,row.names=F,col.names=F,sep="")
	}
}

generate_snp<-function(prefix,i,if_ini=FALSE){
	if(if_ini==TRUE){
		prefixs=paste(prefix,"_",sep="")
	}
	else{
		prefixs=prefix
	}
	if(i==1){
		file.copy(paste(prefix,".snp",sep=""),paste(prefixs,"0.snp",sep=""))
	}
	else{
		file.copy(paste(prefix,".snp",sep=""),paste(prefixs,"1.snp",sep=""))
	}
}

get_results<-function(prefix,bin_dir,if_ini=FALSE){
	fst_evec<-generate_evec(prefix,bin_dir)
	data_matrix<-read_evec(fst_evec)
	if(is.na(get_p(paste(prefix,".tw",sep="")))==FALSE&&get_p(paste(prefix,".tw",sep=""))<=0.01){
		sep_file(prefix,data_matrix,if_ini)
		if(if_ini==TRUE){
			prefix<-paste(prefix,"_",sep="")
		}
		else{
			file.remove(paste(prefix,c(".tw",".ind",".geno",".snp",".pop",".evec",".eval",".log"),sep=""))
		}
		
		left_prefix<-paste(prefix,"0",sep="")
		right_prefix<-paste(prefix,"1",sep="")
		get_results(left_prefix,bin_dir)
		get_results(right_prefix,bin_dir)
	}
}

bin_to_dec<-function(prefix,dir,outfile=NULL){
	dir<-gsub('\\/$','',dir)
	prefix<-gsub('.*\\/','',prefix)
	file_list<-list.files(dir)
#file_list<-file_list[grepl(prefix,file_list)]
#file_list<-file_list[grepl(".ind",file_list)]
	file_list<-file_list[grep(paste0(prefix,'.*[0-1]\\.ind'),file_list)]
#bin_num<-unlist(strsplit(unlist(strsplit(file_list,prefix)),".ind"))
	bin<-gsub(paste0(prefix,"_(.*)\\.ind"),"\\1",file_list)
#bin<-strsplit(bin_num,"")
	n<-length(bin)
	nlen<-1:n
	n_num<-1:n
	for(i in 1:n){
		bin_tmp<-as.numeric(bin[i])
		nlen[i]<-length(bin_tmp)
	}
	nmax<-max(nlen)
	dir<-gsub('(\\/)?$','/',dir)
	inds<-read.table(paste(dir,file_list[1],sep=""))[[1]]
	inds<-as.character(inds)
	cluster<-rep(n_num[1],length(inds))
	for(i in 2:n){
		ind_tmp<-read.table(paste(dir,file_list[i],sep=""))[[1]]
		ind_tmp<-as.character(ind_tmp)
		clu_tmp<-rep(n_num[i],length(ind_tmp))
		inds<-c(inds,ind_tmp)
		cluster<-c(cluster,clu_tmp)
	}
	index<-unlist(strsplit(inds,"SAMPLE"))
	index<-as.numeric(index[which(index!="")])
	dt<-data.frame(inds,cluster,index)
	order_dt<-dt[order(index),][,c(1,2)]
	if(is.null(outfile)){
		write.table(order_dt,paste(dir,prefix,"_res.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
	else{
		write.table(order_dt,outfile,row.names=FALSE,col.names=FALSE,quote=FALSE)
	}
}

########################################## execute part
########################################## 
execute<-function(indir,bin_dir,outfile=NULL){
	if(!file.info(indir)[['isdir']])stop('Not a valid directory!')
	file_list<-list.files(indir)
	if(length(file_list)!=3) stop('File number should be 3!')
	if(!any(grepl('\\.ind',file_list))&&!any(grepl('\\.geno',file_list))&&!any(grepl('\\.snp',file_list))) stop('Lack of at least one of the following files: 1.ind, 2.geno, 3.snp')
	indir<-gsub('(\\/)?$','',indir)
	prefix<-paste0(indir,'/',gsub('\\.ind','',file_list[grep('\\.ind',file_list)][1]))


	get_results(prefix,bin_dir,TRUE)

	if(is.null(outfile)){
		bin_to_dec(prefix,indir)
	}
	else{
		bin_to_dec(prefix,indir,outfile)
	}
}
args<-commandArgs(T)
indir<-args[1]
bin_dir<-args[2]
if(length(args)==3){
	execute(indir,bin_dir,args[3])
}else{
	execute(indir,bin_dir)
}
q(save='no')
