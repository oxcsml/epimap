mkdir website/$1

python3 reinflate.py \
  `ls fits/Rmap-$1*_Rt.csv` `ls fits/Rmap-$1*_Cproj.csv` \
  website/$1/Rt.csv website/$1/Cproj.csv
