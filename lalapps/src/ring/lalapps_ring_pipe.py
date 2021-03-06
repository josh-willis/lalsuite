"""
ring_pipe.py - standalone ring pipeline driver script

This script produces the condor submit and dag files to run
the standalone ring code on LIGO data
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

# import standard modules
import sys, os
import getopt, re, string
import tempfile
import ConfigParser

# import the modules we need to build the pipeline
from glue import pipeline
import ring

def usage():
  msg = """\
Usage: lalapps_ring_pipe [options]

  -h, --help               display this message
  -v, --version            print version information and exit
  -u, --user-tag TAG       tag the job with TAG (overrides value in ini file)

  -d, --datafind           run LSCdataFind to create frame cache files
  -r, --ring              run lalapps_inspiral on the first IFO

  -j, --injections FILE    add simulated bursts from FILE

  -p, --playground-only    only create chunks that overlap with playground
  -P, --priority PRIO      run jobs with condor priority PRIO

  -f, --config-file FILE   use configuration file FILE
  -l, --log-path PATH      directory to write condor log file
"""
  print >> sys.stderr, msg

# pasrse the command line options to figure out what we should do
shortop = "hvdrj:u:P:f:l:"
longop = [
  "help",
  "version",
  "datafind",
  "ring",
  "injections=",
  "playground-only",
  "user-tag=",
  "priority=",
  "config-file=",
  "log-path="
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

config_file = None
do_datafind = None
do_ring = None
inj_file = None
usertag = None
playground_only = 0
condor_prio = None
config_file = None
log_path = None

for o, a in opts:
  if o in ("-v", "--version"):
    print "$Id$"
    sys.exit(0)
  elif o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-d", "--datafind"):
    do_datafind = 1
  elif o in ("-r", "--ring"):
    do_ring = 1
  elif o in ("-j", "--injections"):
    inj_file = a
  elif o in ("-u", "--user-tag"):
    usertag = a
  elif o in ("-p", "--playground-only"):
    playground_only = 1
  elif o in ("-P", "--priority"):
    condor_prio = a
  elif o in ("-f", "--config-file"):
    config_file = a
  elif o in ("-l", "--log-path"):
    log_path = a
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if not config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use --config-file FILE to specify location."
  sys.exit(1)

if not log_path:
  print >> sys.stderr, "No log file path specified."
  print >> sys.stderr, "Use --log-path PATH to specify a location."
  sys.exit(1)

# try and make a directory to store the cache files and job logs
try: os.mkdir('cache')
except: pass
try: os.mkdir('logs')
except: pass

# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(config_file)

# if a usertag has been specified, override the config file
if usertag:
  cp.set('pipeline','user-tag',usertag)
else:
  try:
    usertag = string.strip(cp.get('pipeline','user-tag'))
  except:
    usertag = None

# create a log file that the Condor jobs will write to
basename = re.sub(r'\.ini',r'',config_file)
tempfile.tempdir = log_path
if usertag:
  tempfile.template = basename + '.' + usertag + '.dag.log.'
else:
  tempfile.template = basename + '.dag.log.'
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile)
if usertag:
  dag.set_dag_file(basename + '.' + usertag + '.dag')
else:
  dag.set_dag_file(basename + '.dag')

# create the Condor jobs that will be used in the DAG
df_job = pipeline.LSCDataFindJob('cache','logs',cp)
ring_job = ring.RingJob(cp)

# set better submit file names than the default
if usertag:
  subsuffix = '.' + usertag + '.sub'
else:
  subsuffix = '.sub'
df_job.set_sub_file( basename + '.datafind'+ subsuffix )
ring_job.set_sub_file( basename + '.ring' + subsuffix )

# set the usertag in the jobs
if usertag:
  ring_job.add_opt('user-tag',usertag)

# add the injections
if inj_file:
  ring_job.add_opt('injection-file',inj_file)

# set the condor job priority
if condor_prio:
  df_job.add_condor_cmd('priority',condor_prio)
  ring_job.add_condor_cmd('priority',condor_prio)

# get the pad and chunk lengths from the values in the ini file
pad = 0
length = int(cp.get('data','chunk-length'))
n = int(cp.get('ring', 'segment-duration'))
overlap = n / 2

# read science segs that are greater or equal to a chunk from the input file
data = pipeline.ScienceData()
data.read(cp.get('input','segments'),length)

# create the chunks from the science segments
data.make_chunks(length,overlap,playground_only,0,overlap/2)
#data.make_chunks_from_unused(
#  length,overlap/2,playground_only,overlap/2,0,overlap/2)

# get the order of the ifos to filter
ifo = cp.get('pipeline','ifo')
chan = cp.get('input', 'channel')
ifo_snr = cp.get('pipeline','threshold')

# create all the LSCdataFind jobs to run in sequence
prev_df = None

for seg in data:
  # find all the data
  df = pipeline.LSCDataFindNode(df_job)
  df.set_start(seg.start() - pad)
  df.set_end(seg.end() + pad)
  df.set_observatory(ifo[0])
  if prev_df: 
    df.add_parent(prev_df)

  if do_datafind:
    dag.add_node(df)

  prev_df = df

  for chunk in seg:
    ring_node = ring.RingNode(ring_job)
    ring_node.set_start(chunk.start())
    ring_node.set_end(chunk.end())
    ring_node.set_ifo(ifo)
    ring_node.set_cache(df.get_output())
    ring_node.add_var_opt('threshold',ifo_snr)
    ring_node.add_var_opt('channel-name',ifo+':'+chan)
    if do_datafind: 
      ring_node.add_parent(df)
    if do_ring: 
      dag.add_node(ring_node)

# write out the DAG
dag.write_sub_files()
dag.write_dag()

# write a message telling the user that the DAG has been written
print "\nCreated a DAG file which can be submitted by executing"
print "\n   condor_submit_dag", dag.get_dag_file()
print """\nfrom a condor submit machine (e.g. hydra.phys.uwm.edu)\n
If you are running LSCdataFind jobs, do not forget to initialize your grid 
proxy certificate on the condor submit machine by running the commands

  unset X509_USER_PROXY
  grid-proxy-init -hours 72

Enter your pass phrase when promted. The proxy will be valid for 72 hours. 
If you expect the LSCdataFind jobs to take longer to complete, increase the
time specified in the -hours option to grid-proxy-init. You can check that 
the grid proxy has been sucessfully created by executing the command:

  grid-cert-info -all -file /tmp/x509up_u`id -u`

This will also give the expiry time of the proxy. You should also make sure
that the environment variable LSC_DATAFIND_SERVER is set the hostname and
optional port of server to query. For example on the UWM medusa cluster this
you should use

  export LSC_DATAFIND_SERVER=dataserver.phys.uwm.edu

Contact the administrator of your cluster to find the hostname and port of the
LSCdataFind server.
"""

# write out a log file for this script
if usertag:
  log_fh = open(basename + '.pipeline.' + usertag + '.log', 'w')
else:
  log_fh = open(basename + '.pipeline.log', 'w')
  
# FIXME: the following code uses obsolete CVS ID tags.
# It should be modified to use git version information.
log_fh.write( "$Id$" + "\n\n" )
log_fh.write( "Invoked with arguments:\n" )
for o, a in opts:
  log_fh.write( o + ' ' + a + '\n' )
log_fh.write( "\n" )
log_fh.write( "Parsed " + str(len(data)) + " science segments\n" )
total_data = 0
for seg in data:
  for chunk in seg:
    total_data += len(chunk)
print >> log_fh, "total data =", total_data

print >> log_fh, "\n===========================================\n"
print >> log_fh, data
for seg in data:
  print >> log_fh, seg
  for chunk in seg:
    print >> log_fh, chunk

sys.exit(0)

