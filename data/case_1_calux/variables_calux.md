_This file describes the variable to the data file:_ 

## Data file "zonmw_calux_all_data.csv"

# Variables:

## "Study"	
This variable discribes to which study the data belongs to. Two seprate studies were performed: 							
 - Study 1: nonCoded chemicals were used, in 
 - Study 2: Chemicals were coded								

Type of data: FACTOR

## "Laboratory"	
The tests were performed in Two laboratories: BDS: BioDetection B.V. Amsterdam and HU: Hogeschool Utrecht, lectoraat Innovative Testing
Type of data: DICHOTOME, FACTOR/2

## "Class"	
This indicates the Genotoxic Class to which the chemical belongs, based on IN vivo classification:  									
 - Class = 1: True Positives
 - Class = 2: False Positives 
 - Class = 3: True Negatives

Type data: FACTOR/3: 1, 2 or 3

# "Name"	
The generic name/chemical name of the compound used.																		
Type of data: FACTOR/115: 1 - 115 (or)

# "Nr"
The ID  number of the chemical		
FACTOR/115: 1 - 115 (or)

# "p53.Ifmax"	
The maximum induced fluorescence observed 																	
Type of data: NUMERIC

# "p53.cyto.MEC"																	
Type of data: NUMERIC

# "S9 mix"	
S9 mix was added to some tests: 0=no S9 mix, 1=S9 mix added																	
Type of data: DICHOTOME, FACTOR/2: 1 or 2																		
