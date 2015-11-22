#!/bin/bash

#echo -e string desc=\"$(grep -o -h 'Options::get(.*)' *pp | grep -v shortopt | cut -d '(' -f2 | cut -d ',' -f 1-2 | sort  | uniq | tr -d '0,' | sed "s/'[a-zA-Z]/-&/g" | sed "s/\"[a-zA-Z]/--&/g" | tr -d "\"'" | sed "s/^ //g"| tr '\n' "\\n")\" 

#echo string desc=$(grep -o -h 'Options::get(.*)' *pp | grep -v shortopt | cut -d '(' -f2 | cut -d ',' -f 1-2 | sort  | uniq | tr -d '0,' | sed "s/'[a-zA-Z]/-&/g" | sed "s/\"[a-zA-Z]/--&/g" | tr -d "\"'" | sed "s/^ //g" ) \;
grep -o -h 'Options::get(.*)' *pp | grep -v shortopt | cut -d '(' -f2 | cut -d ',' -f 1-2 | sort  | uniq | tr -d '0,' | sed "s/'[a-zA-Z]/-&/g" | sed "s/\"[a-zA-Z]/--&/g" | tr -d "\"'" | sed "s/^ //g" 
