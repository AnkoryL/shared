#nohup Rscript /data/ivan/rscript/HCMP_miokard/star_rsem_hcmp_hcmp_and_control_geo.r > /data/ivan/nohup_logs/star_rsem_hcmp_hcmp_and_control_geo.r.log 2>&1 &
source("/data/ivan/rscript/a_few_common_functions.r")
suppressMessages(library(ShortRead))

path_input<-"/data/ivan/RNA_seq/HCMP_miokard_GEO/raw"

main_output_folder<-"/data/ivan/RNA_seq/HCMP_miokard_GEO"

output_folder_trim<-paste0(main_output_folder,"/","fastq_trimmed")
output_folder_alignments<-paste0(main_output_folder,"/","aligned")

dir.create(main_output_folder, showWarnings = FALSE)
# dir.create(input_folder, showWarnings = FALSE)
dir.create(output_folder_alignments, showWarnings = FALSE)
dir.create(output_folder_trim, showWarnings = FALSE)
tempfile_adress<-"/data/ivan/RNA_seq/HCMP_miokard_GEO/temp_zcat"

# this needs to remain as a string at all times
num_of_lines_for_fastq_mean_and_sd<-"100000"

qmax<-75
adapter_removal<-"AdapterRemoval"
rsem<-c("/data/ivan/tools/RSEM/rsem-calculate-expression")
star<-c("/data/ivan/tools/STAR-2.7.10b/bin/Linux_x86_64")
index_file<-c("/data/ivan/references/human/star_rsem_indexes/rsem_index")
output_file_fpkm<-c("/data/ivan/RNA_seq/HCMP_miokard_GEO/hcmp_and_control_all_sets_hcmp_geo_rsem_fpkm.txt")
output_file_counts<-c("/data/ivan/RNA_seq/HCMP_miokard_GEO/hcmp_and_control_all_sets_hcmp_geo_rsem_counts.txt")


input_list<-list.files(path_input)
patterns_pair<- c("_1","_2")
pattern_input<-"\\.fastq\\.gz"
pattern_output<-".fastq.gz"
filtered_list<-grep(pattern_input,input_list,value=TRUE)
input_name<-gsub(pattern_input,"",filtered_list)
input_name<-gsub(patterns_pair[1],"",input_name)
input_name<-gsub(patterns_pair[2],"",input_name)
unique_input_name<-unique(input_name)
# dry<-TRUE
dry<-FALSE
debug_mode<-TRUE
do_test<-TRUE
# do_only_one<-TRUE
do_only_one<-FALSE

# trimming
trimming<-TRUE
pattern_trim<-"\\.truncated"
pattern_fastq<-"\\.fastq"

output_trimmed<-list.files(output_folder_trim)
filtered_output_trimmed<-grep(pattern_trim,output_trimmed,value=TRUE)
executed<-0

unique_id_trimmed_singletons<-grep(".singleton",filtered_output_trimmed)
unique_id_trimmed<-filtered_output_trimmed[unique_id_trimmed_singletons*-1]

unique_id_trimmed<-gsub(".gz","",unique_id_trimmed)

unique_id_trimmed<-gsub(".pair1","",unique_id_trimmed)
unique_id_trimmed<-gsub(".pair2","",unique_id_trimmed)
unique_id_trimmed<-gsub(".fastq","",unique_id_trimmed)
unique_id_trimmed<-gsub(".truncated","",unique_id_trimmed)

for (i in 1:length(unique_input_name)) {
  #identify whether input file is paired or not by mathcing unique part of the name to the name vector 
  match_name<-input_name==unique_input_name[i]
  how_many_matches<-sum(match_name)
  if (how_many_matches==1) {
    paired<-FALSE	
  } else if(how_many_matches==2) {
    paired<-TRUE	
  } else {
    print(paste0("This shouldn't have happened - ",how_many_matches," matches found"))
    print(unique_input_name[i]) 
    print(paste0(input_name[match_name],sep="",collapse=" "))
    break
  }
  if (unique_input_name[i] %in% unique_id_trimmed) {next}
  #/identify
  
  # inner test - both are present in the original list
  if (paired) {
    target_file<-paste0(unique_input_name[i],patterns_pair,pattern_output)
    target_file_paired_out<-paste0(unique_input_name[i],pattern_output)
    # print(target_file)
  } else {
    target_file<-paste0(unique_input_name[i],pattern_output)
    target_file_paired_out<-paste0(unique_input_name[i],pattern_output)
    # print(target_file)
  } 
  if (do_test) {
    if (sum(target_file %in% filtered_list)==length(target_file)) {
      # print(paste0("Seems ok, constructed filename is present in the original list ",i,"/",length(unique_input_name)))
    } else {
      print(paste0("Error, constructed filename is NOT present in the original list",i,"/",length(unique_input_name)))
      print(target_file)
    }
  }
  
  #/inner test
  
  #creating and running the command
  
  if (!paired) {
    com<-paste(adapter_removal," --file1 ",path_input,"/",target_file[1]," --basename ",output_folder_trim,"/",target_file_paired_out, " --gzip --threads 32 --trimns --trimqualities --qualitymax ",qmax, sep="",collapse=NULL)
    # AdapterRemoval --file1 /data/rats_transcriptomes/raw/170322_HSGA.Slominsky_RNA_2017.Slo_G6_II.fastq --basename output_single --trimns --trimqualities
  } else {
    com<-paste(adapter_removal," --file1 ",path_input,"/",target_file[1]," --file2 ",path_input,"/",target_file[2]," --basename ",output_folder_trim,"/",target_file_paired_out, " --gzip --threads 32 --trimns --trimqualities --qualitymax ",qmax, sep="",collapse=NULL)
    # AdapterRemoval --file1 /data/rats_transcriptomes/raw/170322_HSGA.Slominsky_RNA_2017.Slo_G6_II.fastq --basename output_single --trimns --trimqualities
  }		
  
  if (do_only_one & (executed>1)) {break}
  if(debug_mode) {
    print(com)
  }
  if (!dry) {
    system(com)
    
    executed<-executed+1
  }
}
# done with trimming


# rename files to make sure the patterns of names are correct and uniform
# output_trimmed_suffixes<-gsub("SRR[0-9]+\\.","",output_trimmed)
# unique(output_trimmed_suffixes)

output_trimmed<-list.files(output_folder_trim)


remove_all_but_last_unique_instance<-function(input,pattern="\\.") {
  input_splat<-unlist(strsplit(input,pattern),use.names=FALSE)
  u_test<-unique(input_splat)
  out_numbers<-c()
  for (i in 1:length(u_test)) {
    # print(u_test[i])
    where<-grep(u_test[i],input_splat)
    where_last_encounter<-max(where)
    out_numbers<-c(out_numbers,where_last_encounter)
  }
  out_numbers<-out_numbers[order(out_numbers)]
  output<-input_splat[out_numbers]
  
  output<-paste(output,collapse=".",sep=".")
  return(output)
}
output_trimmed_fixed<-unname(sapply(output_trimmed,remove_all_but_last_unique_instance,pattern="\\."))

output_trimmed_fixed_full_path<-paste0(output_folder_trim,"/",output_trimmed_fixed)
output_trimmed_full_path<-paste0(output_folder_trim,"/",output_trimmed)


for (i in 1:length(output_trimmed_full_path)) {
  if (output_trimmed_full_path[i]==output_trimmed_fixed_full_path[i]) {
    next
  } else {
    if (debug_mode) {
      print(paste0("Renaming:",output_trimmed_full_path[i]," to ",output_trimmed_fixed_full_path[i]))
    }
    if (!dry) {
      file.rename(output_trimmed_full_path[i],output_trimmed_fixed_full_path[i])
    }
  }
}

# output_trimmed_suffixes<-gsub("SRR[0-9]+\\.","",output_trimmed_fixed)
# unique(output_trimmed_suffixes)
#/rename


# alignment
pattern_aligned<-"\\.genes\\.results"
output_aligned<-list.files(output_folder_alignments)

filtered_output_aligned<-grep(pattern_aligned,output_aligned,value=TRUE)
# filtered_output_aligned<-grep(".fastq",filtered_output_aligned,value=TRUE)
# filtered_output_aligned<-grep(".truncated",filtered_output_aligned,value=TRUE)


id_aligned<-gsub(pattern_aligned,"",filtered_output_aligned)
id_aligned<-gsub(pattern_trim,"",id_aligned)

id_aligned<-gsub(pattern_fastq,"",id_aligned)

unique_id_aligned<-unique(id_aligned)

patterns_pair_trimmed<-c(".pair1",".pair2")
pattern_trim<-".truncated"
pattern_gzip<-".gz"
pattern_aligned<-".bam"

input_list<-list.files(output_folder_trim)
input_list_fixed<-unname(sapply(input_list,remove_all_but_last_unique_instance,pattern="\\.fastq."))
input_list_fixed_full_path<-paste0(output_folder_trim,"/",input_list_fixed)
input_list_full_path<-paste0(output_folder_trim,"/",input_list)
for (i in 1:length(input_list_full_path)) {
  if (input_list_full_path[i]==input_list_fixed_full_path[i]) {
    next
  } else {
    if (debug_mode) {
      print(paste0("Renaming:",input_list_full_path[i]," to ",input_list_fixed_full_path[i]))
    }
    if (!dry) {
      file.rename(input_list_full_path[i],input_list_fixed_full_path[i])
    }
  }
}

input_list<-list.files(output_folder_trim)

filtered_list<-grep(pattern_trim,input_list,value=TRUE)
if (paired) {
  singleton_list<-grep("singleton",filtered_list) 
  filtered_list<-filtered_list[-1*singleton_list]
}
input_name<-gsub(pattern_trim,"",filtered_list)
input_name<-gsub(patterns_pair_trimmed[1],"",input_name)
input_name<-gsub(patterns_pair_trimmed[2],"",input_name)
input_name<-gsub(pattern_gzip,"",input_name)
input_name<-gsub(pattern_fastq,"",input_name)

unique_input_name<-unique(input_name)

executed<-0

for (i in 1:length(unique_input_name)) {
  if (unique_input_name[i] %in% unique_id_aligned) {next}
  
  #identify whether input file is paired or not by mathcing unique part of the name to the name vector 
  
  match_name<-input_name==unique_input_name[i]
  how_many_matches<-sum(match_name)
  
  if (how_many_matches==1) {
    paired<-FALSE	
  } else if (how_many_matches==2) {
    paired<-TRUE	
  } else {
    print(paste0("This shouldn't have happened - ",how_many_matches," matches found"))
    print(unique_input_name[i]) 
    print(paste0(input_name[match_name],sep="",collapse=" "))
    break
  }
  #/identify
  
  
  # inner test - both are present in the original list
  if (paired) {
    target_file<-paste0(unique_input_name[i],patterns_pair_trimmed,pattern_trim,pattern_gzip)
    target_file_path<-paste0(output_folder_trim,"/",target_file)
    target_file_paired_out<-paste0(unique_input_name[i])
    alignments_fullpath<-paste0(output_folder_alignments,"/",target_file_paired_out)
    # print(target_file)
  } else {
    target_file<-paste0(unique_input_name[i],pattern_trim,pattern_gzip)
    target_file_path<-paste0(output_folder_trim,"/",target_file)
    target_file_paired_out<-paste0(unique_input_name[i])
    alignments_fullpath<-paste0(output_folder_alignments,"/",target_file_paired_out)
    # print(target_file)
  } 
  if (do_test) {
    if (sum(target_file %in% filtered_list)==length(target_file)) {
      # print(paste0("Seems ok, constructed filename is present in the original list ",i,"/",length(unique_input_name)))
    } else {
      print(paste0("Error, constructed filename is NOT present in the original list ",i,"/",length(unique_input_name)))
      print(target_file)
    }
  }
  # getting average read length and standart deviation
  if (dry) {
    reads_length_mean<-100
    reads_length_sd<-20
  } else {
    
    com_zcat<-paste0("zcat ",target_file_path[1]," | head -",num_of_lines_for_fastq_mean_and_sd," > ",tempfile_adress)
    
    if (debug_mode) {
      print(com_zcat)
    }
    
    system(com_zcat)
    
    fastq_file<-readFastq(tempfile_adress)
    reads_length_mean<-mean(width(sread(fastq_file)))
    reads_length_sd<-sd(width(sread(fastq_file)))
    fastq_file<-c()
    unlink(tempfile_adress)
  }
  if (paired) {
    
    # print("check")
    com<-paste0(rsem," --output-genome-bam --star-gzipped-read-file --star --star-path ",
                star, " \\\n ", "--phred33-quals"," \\\n ","--fragment-length-mean ",
                reads_length_mean," \\\n ", " --fragment-length-sd ",reads_length_sd,
                " \\\n ", "-p 64 --calc-ci --ci-memory 61024 "," --paired-end ", 
                target_file_path[1]," \\\n ",target_file_path[2]," \\\n ",index_file," \\\n " ,alignments_fullpath)
  } else {
    # print("check")
    
    com<-paste0(rsem," --output-genome-bam --star-gzipped-read-file --star --star-path ",
                star, " \\\n ", "--phred33-quals"," \\\n ","--fragment-length-mean ",
                reads_length_mean," \\\n ", " --fragment-length-sd ",reads_length_sd,
                " \\\n ", "-p 64 --calc-ci --ci-memory 61024 ", 
                target_file_path[1]," \\\n ",index_file," \\\n " ,alignments_fullpath)
  }
  
  if (do_only_one & (executed>1)) {break}
  if(debug_mode) {
    print(com)
  }
  if (!dry) {
    system(com)
  }
  executed<-executed+1
}


# input_list_aligned<-list.files(output_folder_alignments)
# filelist_aligned_results<-grep(".genes.results",input_list_aligned,value=TRUE)
# filelist_aligned_results_paths<-paste0(output_folder_alignments,"/",filelist_aligned_results)

# input_frame<-data.frame(Gene_id=character(), Transcript_id=character(), Length=numeric())
# for (i in 1:length(filelist_aligned_results_paths)) {
# rsem_results<-read.table(filelist_aligned_results_paths[i],sep="\t",quote="",header=TRUE)
# if (i==1) {
# input_matrix<-matrix(,nrow=nrow(rsem_results),ncol=length(filelist_aligned_results_paths))
# rownames(input_matrix)<-rsem_results[,1]
# input_frame<-rsem_results[,1:3]
# }
# input_matrix[,i]<-unlist(rsem_results[,5],use.names=FALSE)
# }

# colnames(input_matrix)<-filelist_aligned_results
# if (!dry) {
# write.table(input_matrix,output_file_counts,quote = FALSE,sep = "\t",row.names = TRUE)
# }

# input_frame<-data.frame(Gene_id=character(), Transcript_id=character(), Length=numeric())
# for (i in 1:length(filelist_aligned_results_paths)) {
# rsem_results<-read.table(filelist_aligned_results_paths[i],sep="\t",quote="",header=TRUE)
# if (i==1) {
# input_matrix<-matrix(,nrow=nrow(rsem_results),ncol=length(filelist_aligned_results_paths))
# rownames(input_matrix)<-rsem_results[,1]
# input_frame<-rsem_results[,1:3]
# }
# input_matrix[,i]<-unlist(rsem_results[,7],use.names=FALSE)
# }

# colnames(input_matrix)<-filelist_aligned_results
# if (!dry) {
# write.table(input_matrix,output_file_fpkm,quote = FALSE,sep = "\t",row.names = TRUE)
# }