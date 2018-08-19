library(data.table)
library(foreach)
library(doMC)
registerDoMC(15)
library(stringr)

# Need to run on durga

# Published GWAS:
in_dir = '/users/mgloud/projects/gwas/data/munged/'
out_dir = 'processed_data/count_trait/count_trait/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

in_fn_list = list.files(in_dir,pattern='gz$',full.names=TRUE)
message(length(in_fn_list),' published GWAS.')

count_trait = function(in_fn_list,out_dir){
	n_trait = foreach(fn = in_fn_list,.combine='c')%do%{
		print(fn)
		command = sprintf('zcat %s | cut -f1-1 | uniq',fn)
		result = system(command,intern=TRUE)
		if (str_detect(result[1],'rsid')){
			traits = str_split_fixed(basename(fn),'_',4)[,2]
			n_trait = 1
		} else {
			result = result[result!='trait']
			traits = unique(result)
			n_trait = length(traits)
		}
		out_fn = sprintf('%s/%s.txt',out_dir,basename(fn))
		write.table(traits,out_fn,sep='\t',col.names=FALSE,row.names=FALSE,quote=FALSE)
		rm(result)
		return(n_trait)
	}
	return(n_trait)
}

n_trait = count_trait(in_fn_list,out_dir)

fn_list = list.files(out_dir,full.names=TRUE)

trait = foreach(fn = fn_list,.combine='rbind')%dopar%{
	fread(fn,header=FALSE)
}
trait = unname(unlist(trait))
length(unique(trait)) # 789

# UKBB:
in_dir = '/users/mgloud/projects/gwas/data/ukbb/'
in_fn_list = list.files(in_dir,pattern='txt.gz$',full.names=TRUE)
message(length(in_fn_list)+1,' UKBB GWAS.')

n_ukbb_trait = foreach(fn = in_fn_list,.combine='c')%do%{
	print(fn)
	command = sprintf('zcat %s | cut -f1-1 | uniq',fn)
	result = system(command,intern=TRUE)
	split_result = unlist(str_split(result,pattern='\n'))
	split_result = split_result[split_result!='trait']
	return(split_result)
}
n_ukbb_trait = c("Q28","Q21","Q18","Q82","Q38","K85","K46","K50","K65","K25","K29","K86", "K56","K12","K22","K60","K40","K57","K91","K62","K00","K75","K90","K05",  "K51","K83","K74","K21","K58","K35","K52","K42","K04","K61","K66","K81",  "K06","K92","K08","K59","K37","K09","K31","K26","K70","K11","K44","K10",  "K82","K13","K02","K01","K30","K07","K63","K80","K20","K14","K27","K41",  "K43","K76","K55","E66","E10","E21","E03","E27","E23","E07","E14","E04",  "E87","E16","E83","E22","E11","E86","E05","H80","H34","H16","H01","H43",  "H71","H66","H40","H35","H65","H00","H72","H74","H04","H53","H26","H02",  "H59","H91","H83","H52","H27","H69","H33","H57","H50","H92","H11","H05", "H90","H61","H49","H18","H60","H73","H25","H93","H81","Z01","Z41","Z44", "Z04","Z08","Z87","Z32","Z43","Z13","Z47","Z35","Z00","Z39","Z50","Z76", "Z34","Z45","Z48","Z71","Z12","Z09","Z03","Z40","Z42","Z31","Z36","Z52", "Z85","Z53","Z80","Z30","Z51","Z11","Z46","N81","N49","N45","N42","N88", "N82","N13","N50","N43","N86","N36","N40","N30","N99","N84","N23","N61", "N62","N90","N39","N89","N87","N48","N05","N02","N76","N19","N17","N35", "N73","N20","N28","N31","N71","N34","N63","N95","N70","N64","N41","N92", "N83","N94","N72","N85","N97","N93","N12","N10","N32","N75","N21","N80", "N47","N60","N18","B50","B07","A87","B18","A09","B37","B02","A63","A08", "B34","A04","A41","S72","S12","S62","S39","T83","T40","S64","S56","T78", "S43","T86","S27","S52","T85","T14","S05","S22","T18","T50","S13","S46", "T43","S93","S67","S32","S51","S86","T87","S92","T17","S00","S76","T39", "S61","S01","S91","S63","S06","S09","S81","S83","S82","S20","S30","S69", "S80","T82","S68","S42","T84","T88","T81","S31","S89","S66","S60","T42", "S02","T79","F30","F25","F99","F31","F45","F41","F52","F43","F20","F10", "F23","F32","F60","F33","M00","M24","M77","M50","M13","M95","M20","M47", "M81","M18","M79","M21","M06","M72","M34","M17","M67","M46","M10","M75", "M22","M19","M16","M05","M41","M23","M31","M85","M25","M60","M71","M48", "M62","M86","M65","M35","M94","M15","M76","M84","M70","M87","M53","M80", "M89","M96","M51","M54","M66","M43","M32","M45","D16","D47","D09","D32", "D26","D38","D17","D36","D37","D07","D05","D86","D69","D64","D46","D25", "D75","D30","D12","D23","D34","D41","D21","D33","D18","D68","D06","D61", "D27","D70","D11","D22","D39","D24","D44","D13","D28","D14","D35","D03", "D89","D48","D50","D45","D10","D04","C82","C50","C21","C32","C71","C77", "C64","C85","C02","C69","C81","C91","C90","C80","C79","C73","C15","C43", "C56","C54","C62","C16","C78","C49","C19","C83","C61","C09","C92","C20", "C22","C25","C45","C44","C53","C34","C67","C18","G45","G24","G54","G50", "G35","G82","G56","G20","G40","G47","G93","G57","G51","G81","G58","G43", "G61","G70","G62","G95","G44","G37","O46","O16","O69","O64","O48","O00", "O60","O63","O03","O44","O80","O68","O73","O14","O04","O34","O82","O70", "O32","O41","O72","O99","O13","O35","O90","O62","O26","O66","O24","O23", "O42","O36","O47","O02","O75","O21","O81","O20","J01","J31","J33","J84", "J98","J34","J36","J44","J32","J38","J86","J40","J39","J46","J43","J22", "J02","J45","J35","J06","J13","J15","J90","J93","J18","J47","J96","J03", "L29","L90","L03","L82","L05","L30","L50","L91","L81","L27","L72","L92", "L85","L53","L02","L57","L73","L28","L60","L97","L08","L98","L40","L43", "R40","R87","R68","R93","R09","R56","R73","R00","R53","R51","R15","R54", "R29","R18","R69","R55","R06","R22","R10","R20","R32","R47","R52","R19", "R91","R49","R13","R21","R60","R26","R50","R39","R14","R05","R07","R42", "R61","R90","R79","R86","R12","R59","R35","R17","R31","R30","R25","R23", "R63","R11","R76","R41","R04","R33","R94","I20","I21","I25","I82","I35","I74","I48","I60","I61","I85","I10","I08","I65","I26","I67","I33","I63","I84","I31","I42","I72","I70","I71","I86","I51","I89","I78","I47","I24","I95","I45","I64","I50","I73","I80","I77","I62","I87","I22","I12","I05","I30","I46","I49","I83","I44","I34")
length(unique(n_ukbb_trait)) # 642