sed -e '/^>/s/$/@/' -e 's/^>/#/' $1  |\
 tr -d '\n' |\
 tr "#" "\n" |\
 tr "@" " " |\
 sed '/^$/d' |\
 sort -u -t " " -f -k 2,2 |\
 sed -e 's/^/>/' |\
 tr -s \  '\n'
