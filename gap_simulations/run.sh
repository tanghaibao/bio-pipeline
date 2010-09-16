COUNT_RETENTIONS=1
#python plot_sims.py small.runs.a.txt 1 sim 1

python plot_sims.py runs.a.txt $COUNT_RETENTIONS sim 1
python plot_sims.py runs.a.txt $COUNT_RETENTIONS sim 2
python plot_sims.py runs.a.txt $COUNT_RETENTIONS sim 3
python plot_sims.py runs.a.txt $COUNT_RETENTIONS sim 4
python plot_sims.py runs.b.txt $COUNT_RETENTIONS sim 1
python plot_sims.py runs.b.txt $COUNT_RETENTIONS sim 2
python plot_sims.py runs.b.txt $COUNT_RETENTIONS sim 3
python plot_sims.py runs.b.txt $COUNT_RETENTIONS sim 4
python plot_sims.py runs.c.txt $COUNT_RETENTIONS sim 1
python plot_sims.py runs.c.txt $COUNT_RETENTIONS sim 2
python plot_sims.py runs.c.txt $COUNT_RETENTIONS sim 3
python plot_sims.py runs.c.txt $COUNT_RETENTIONS sim 4
