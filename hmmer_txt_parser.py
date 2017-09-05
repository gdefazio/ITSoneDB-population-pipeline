__author__ = "Bruno Fosso"
__version__ = "1.0"
import os
from string import strip


if os.path.exists("hmm_match_data/align_18S_data.txt"):
    print "18S alignment data exist"
else:
    print "No 18S alignment data"
    exit()

if os.path.exists("hmm_match_data/align_5_8S_data.txt"):
    print "5.8S alignment data exist"
else:
    print "No 5.8S alignment data"
    exit()


print "5.8S ALIGNEMENT EXTRACTION"

os.mkdir("hmm_match_data/5_8S_alignment_file")
f = open("hmm_match_data/align_5_8S_data.txt")
lines = f.readlines()
f.close()

count = 0
stop_file = len(lines) - 15
line2acc = {}
line_num = []
while count <= stop_file:
    if lines[count][0:2] == ">>":
        s = map(strip, lines[count].split(" "))
        acc = s[1]
        print acc
        line2acc[count] = acc
        line_num.append(count)
        count += 1
    else:
        count += 1

line_num.sort()
index = 0
stop = len(line_num)
while index < stop:
    if index != stop - 1:
        inizio = line_num[index]
        fine = line_num[index + 1]
        acc = line2acc[inizio]
        #print inizio,fine
        tmp = open("hmm_match_data/5_8S_alignment_file/" + acc + ".align", "w")
        tmp.write("".join(lines[inizio:fine]))
        tmp.close()
        index += 1
    else:
        inizio = line_num[index]
        fine = stop_file
        acc = line2acc[inizio]
        tmp = open("hmm_match_data/5_8S_alignment_file/" + acc + ".align", "w")
        tmp.write("".join(lines[inizio:fine]))
        tmp.close()
        index += 1

print "18S ALIGNEMENT EXTRACTION"
os.mkdir("hmm_match_data/18S_alignment_file")
f = open("hmm_match_data/align_18S_data.txt")
lines = f.readlines()
f.close()

count = 0
stop_file = len(lines) - 15
line2acc = {}
line_num = []
while count <= stop_file:
    if lines[count][0:2] == ">>":
        s = map(strip, lines[count].split(" "))
        acc = s[1]
        #print acc
        line2acc[count] = acc
        line_num.append(count)
        count += 1
    else:
        count += 1

line_num.sort()
index = 0
stop = len(line_num)
while index < stop:
    if index != stop - 1:
        inizio = line_num[index]
        fine = line_num[index + 1]
        acc = line2acc[inizio]
        #print inizio,fine
        tmp = open("hmm_match_data/18S_alignment_file/" + acc + ".align", "w")
        tmp.write("".join(lines[inizio:fine]))
        tmp.close()
        index += 1
    else:
        inizio = line_num[index]
        fine = stop_file
        acc = line2acc[inizio]
        tmp = open("hmm_match_data/18S_alignment_file/" + acc + ".align", "w")
        tmp.write("".join(lines[inizio:fine]))
        tmp.close()
        index += 1

print "FINISH"