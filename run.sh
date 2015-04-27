

for i in {1..2000}
do
   echo "CASE:" $i
   time java -jar tester.jar -exec ./program -seed $i
done	
