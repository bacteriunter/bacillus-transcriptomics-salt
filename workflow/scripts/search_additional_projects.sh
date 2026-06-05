#!/usr/bin/env bash

mkdir -p metadata/study_selection/searches/results

while IFS= read -r query <&3
do
    safe=$(echo "$query" | tr ' ' '_' | tr '/' '_')
    outfile="metadata/study_selection/searches/results/${safe}.csv"

    echo "Searching: $query"

    esearch -db sra -query "$query" \
    | efetch -format runinfo \
    > "$outfile" || true

done 3< metadata/study_selection/searches/search_queries.txt
