#!/bin/bash

# Start a Docker container
docker pull quay.io/felicityallen/selftarget
docker run -d --name selftarget -p 5001:8006 quay.io/felicityallen/selftarget

# Copy local file to docker container and split
docker cp $1 selftarget:/app/indel_prediction/predictor/forecast_batch.txt

# Exec into container
docker exec -it selftarget bash
cd indel_prediction/predictor/

# Split file (10000 lines each including header)
tail -n +2 forecast_batch.txt | split -l 9999 - forecast_batch_split_
for file in forecast_batch_split_*
do
  head -n 1 forecast_batch.txt > tmp_file
  cat "$file" >> tmp_file
  mv -f tmp_file "$file"
done

# Run FORECasT batch mode
for file in forecast_batch_split_a*
do
  python FORECasT.py ${file} ${file}_guideSeq79
done

# Copy output files from container to local host
#for letter in {a..z}
#do
#  docker cp selftarget:/app/indel_prediction/predictor/forecast_batch_split_a${letter}_guideSeq79_predictedindelsummary.txt ./
#  docker cp selftarget:/app/indel_prediction/predictor/forecast_batch_split_a${letter}_guideSeq79_predictedreads.txt ./
#done
