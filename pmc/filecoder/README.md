# Product-Matrix File Coder
Practical file-based encoder/collector/repair programs for product-matrix codes.

Parameters:
* n: total number of nodes
* k: any (k) nodes sufficient to recover all data
* d: any (d) nodes sufficient to repair failed node
* w: size of Galois field in bits

## Running locally
All examples use n=5, k=3, d=4, w=8.
Places encoded chunks in coded_files dir, and helper chunks in helper_files dir.

### Generating Sample Files
10Mb random file:
```
dd if=/dev/urandom of=f10M.dat bs=1024 count=10000
```

### Encoding

```
$ ./encoder_pmc 5 3 4 8 MSR video.mp4
Bytes read      : 41393736
Bytes written   : 68989600
Coding Time [s] : 0.260
IO Time [s]     : 0.140

$ rm video.mp4 && ls coded_files
coded_files/video.mp4_metadata
coded_files/video.mp4_node0_sym0
coded_files/video.mp4_node0_sym1
coded_files/video.mp4_node1_sym0
coded_files/video.mp4_node1_sym1
coded_files/video.mp4_node2_sym0
coded_files/video.mp4_node2_sym1
coded_files/video.mp4_node3_sym0
coded_files/video.mp4_node3_sym1
coded_files/video.mp4_node4_sym0
coded_files/video.mp4_node4_sym1
```

### Node Repair
Generate helper chunks for node 0 failing, from nodes 1, 2, 3, 4:

```
$ rm coded_files/video.mp4_node0*
$ ./repair_pmc help 0 coded_files/video.mp4_metadata 1 2 3 4
n: 5, k: 3, d: 4, w: 8, alpha: 2, B: 6, MBR: 0
msg_size: 41393760
helper nodes: 1, 2, 3, 4, 
*** Generating repair chunk from node 1 to failed node 0...
reading chunk: coded_files/video.mp4_node1_sym0
reading chunk: coded_files/video.mp4_node1_sym1
writing helper chunk: helper_files/video.mp4_helper_chunk_from1_to0
*** Generating repair chunk from node 2 to failed node 0...
reading chunk: coded_files/video.mp4_node2_sym0
reading chunk: coded_files/video.mp4_node2_sym1
writing helper chunk: helper_files/video.mp4_helper_chunk_from2_to0
*** Generating repair chunk from node 3 to failed node 0...
reading chunk: coded_files/video.mp4_node3_sym0
reading chunk: coded_files/video.mp4_node3_sym1
writing helper chunk: helper_files/video.mp4_helper_chunk_from3_to0
*** Generating repair chunk from node 4 to failed node 0...
reading chunk: coded_files/video.mp4_node4_sym0
reading chunk: coded_files/video.mp4_node4_sym1
writing helper chunk: helper_files/video.mp4_helper_chunk_from4_to0
Bytes read      : 55191680
Bytes written   : 27595840
Coding Time [s] : 0.020
IO Time [s]     : 0.080
```

Now use helper chunks to repair node 0:

```
$ ./repair_pmc fix 0 coded_files/video.mp4_metadata 1 2 3 4
n: 5, k: 3, d: 4, w: 8, alpha: 2, B: 6, MBR: 0
msg_size: 41393760
helper nodes: 1, 2, 3, 4, 
*** Repairing failed node 0...
reading chunk: helper_files/video.mp4_helper_chunk_from1_to0
reading chunk: helper_files/video.mp4_helper_chunk_from2_to0
reading chunk: helper_files/video.mp4_helper_chunk_from3_to0
reading chunk: helper_files/video.mp4_helper_chunk_from4_to0
writing repaired chunk: coded_files/video.mp4_node0_sym0
writing repaired chunk: coded_files/video.mp4_node0_sym1
Bytes read      : 27595840
Bytes written   : 13797920
Coding Time [s] : 0.080
IO Time [s]     : 0.040

$ ls coded_files/video.mp4_node0*
coded_files/video.mp4_node0_sym0
coded_files/video.mp4_node0_sym1
```

### Data Collection

Collect from nodes 2, 3, 4 to recover full data:

```
$ ./collector_pmc coded_files/video.mp4_metadata 2 3 4
n: 5, k: 3, d: 4, w: 8, alpha: 2, B: 6, MBR: 0
file size: 41393736, msg_size: 41393760
participating nodes: 2, 3, 4, 
reading chunk: coded_files/video.mp4_node2_sym0
reading chunk: coded_files/video.mp4_node2_sym1
reading chunk: coded_files/video.mp4_node3_sym0
reading chunk: coded_files/video.mp4_node3_sym1
reading chunk: coded_files/video.mp4_node4_sym0
reading chunk: coded_files/video.mp4_node4_sym1
*** Writing decoded file: video.mp4
Bytes read      : 41393760
Bytes written   : 41393736
Coding Time [s] : 0.240
IO Time [s]     : 0.100
```


## Running networked

All functional nodes in the network are "servers". Client operations are:

1. Sending out encoded files for servers to store.
2. Requesting data collection.
3. Requesting node repair.

For (2) and (3), the client recieves (processed) chunks from participating nodes, and performs local operations to recover/repair data.

Clients need the node-metadata file, containing the list of server hosts:ports (details in [doc](../doc/metadata_formats.txt)).

Server required files:
- repair_pmc
- pmc_server.py
- common.py

Encoder client required files:
- encoder_pmc
- pmc_client.py
- common.py

Data-collector client required files:
- collector_pmc
- pmc_client.py
- common.py

Failed node client required files:
- repair_pmc
- pmc_client.py
- common.py


For each server, run: ```python pmc_server.py [port] [node id] [host ip]```

For client operations, see options of ```pmc_client.py```



#### Special instructions for running on Windows:

- Install Cygwin (with gcc & make) to compile.
- Add cygwin/bin to PATH
