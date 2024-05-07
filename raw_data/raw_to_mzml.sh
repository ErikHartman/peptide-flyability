for f in chymotrypsin/*.raw
do
 echo "Processing $f" 
 path=/Users/erikhartman/dev/peptide-flyability/raw_data/${f//".raw"/".mzml"}
 echo "Output in $path"
 mono ../../ThermoRawFileParser1.4.3/ThermoRawFileParser.exe -i=$f -b=$path -f=2
done


for f in trypsin/*.raw
do
 echo "Processing $f" 
 path=/Users/erikhartman/dev/peptide-flyability/raw_data/${f//".raw"/".mzml"}
 echo "Output in $path"
 mono ../../ThermoRawFileParser1.4.3/ThermoRawFileParser.exe -i=$f -b=$path -f=2
done


for f in lysyl_pepsin/*.raw
do
 echo "Processing $f"
 path=/Users/erikhartman/dev/peptide-flyability/raw_data/${f//".raw"/".mzml"}
 echo "Output in $path"
 mono ../../ThermoRawFileParser1.4.3/ThermoRawFileParser.exe -i=$f -b=$path -f=2
done


for f in elastase/*.raw
do
 echo "Processing $f" 
 path=/Users/erikhartman/dev/peptide-flyability/raw_data/${f//".raw"/".mzml"}
 echo "Output in $path"
 mono ../../ThermoRawFileParser1.4.3/ThermoRawFileParser.exe -i=$f -b=$path -f=2
done