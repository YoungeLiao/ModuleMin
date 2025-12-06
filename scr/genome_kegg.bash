### ================= v2 ====================
awk -F'\t' '
BEGIN {
    total=0
    valid=0
}
NR==1{
    
    for(i=1;i<=NF;i++){
        if($i=="Gene ID"||$i=="Genome"||$i=="KO"||$i=="KEGG_name"||$i=="KO_description"||$i=="Level3"){
            cols[++c]=i
            # 打印表头
            printf "%s%s", $i, (c==6?ORS:OFS)
        }
    }
    next
}
{
    total++
    
    is_valid=1
    for(j=1;j<=c;j++){
        if($cols[j]=="-"){
            is_valid=0
            break
        }
    }
    if(is_valid){
        valid++
        # 打印有效行
        for(j=1;j<=c;j++){
            printf "%s%s", $cols[j], (j==c?ORS:OFS)
        }
    }
}
END {
    # 打印统计信息到stderr
    print "Total genes: " total > "/dev/stderr"
    print "Valid genes (no missing annotations): " valid > "/dev/stderr"
    print "Percentage of valid genes: " (valid/total*100) "%" > "/dev/stderr"
}' anno_summary.txt > extracted_data.tsv


# Total genes: 110153
# Valid genes (no missing annotations): 48850
# Percentage of valid genes: 44.3474%


Rscript -e "
library(writexl)
df <- read.table('extracted_data.tsv', sep='\t', header=TRUE, quote='', comment.char='', stringsAsFactors=FALSE)
write_xlsx(df, 'extracted_data.xlsx')
unlink('extracted_data.tsv')
"

### ================= v1 ====================
awk -F'\t' '
NR==1{
    # 找到目标列的位置
    for(i=1;i<=NF;i++){
        if($i=="Gene ID"||$i=="Genome"||$i=="KO"||$i=="KEGG_name"||$i=="KO_description"){
            cols[++c]=i
            # 打印表头
            printf "%s%s", $i, (c==5?ORS:OFS)
        }
    }
    next
}
{
    # 打印数据行
    for(j=1;j<=c;j++){
        printf "%s%s", $cols[j], (j==c?ORS:OFS)
    }
}' anno_summary.txt > extracted_data.tsv

# 转换为xlsx
Rscript -e "
library(writexl)
df <- read.table('extracted_data.tsv', sep='\t', header=TRUE, quote='', comment.char='', stringsAsFactors=FALSE)
write_xlsx(df, 'extracted_data.xlsx')
unlink('extracted_data.tsv')
"
