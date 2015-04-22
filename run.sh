

for i in {1..50}
do
   echo "CASE:" $i
   java -jar tester.jar -exec ./program -seed $i
done	
