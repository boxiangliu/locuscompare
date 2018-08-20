library(data.table)

in_fn = '/srv/persistent/bliu2/rpe/processed_data/finemap/manhattan/2018-06-23_21-47-08_rpe_amd_gtex/Fritsche_sorted_txt_gz_finemap_clpp_status.txt'

read_colocalization = function(fn){
    x = fread(fn)
    return(x)
}

create_table = function(){
    
}

upload_table = function(){
    
}


colocalization = read_colocalization(in_fn)
create_table(table_name)
upload_table(table_name)

