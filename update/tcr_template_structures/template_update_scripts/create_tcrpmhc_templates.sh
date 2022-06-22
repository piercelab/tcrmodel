
#read class1 complex structure and cap them to use as template in TCRmodel
while read f ; do python2.7 /www/tcrmodel/update/tcr_template_structures/template_update_scripts/tcrpmhc1_cap.py ${f}; done < /www/tcrmodel/update/tcr_template_structures/files/class1_list.txt

#read class2 complex structure and cap them to use as template in TCRmodel
while read f ; do python2.7 /www/tcrmodel/update/tcr_template_structures/template_update_scripts/tcrpmhc2_cap.py ${f}; done < /www/tcrmodel/update/tcr_template_structures/files/class2_list.txt
