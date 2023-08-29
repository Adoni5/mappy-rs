//! A multithreaded minimap2 mappy clone in rust.
//! Serves as a drop in for mappy in single threaded mode.
//! Designed for use with readfish https://github.com/LooseLab/readfish/
#![deny(missing_docs)]
#![deny(clippy::missing_docs_in_private_items)]

use crossbeam::channel::{bounded, Receiver, RecvError, Sender};
use crossbeam::queue::ArrayQueue;
use fnv::FnvHashMap;
use itertools::all;
use pyo3::exceptions::{
    PyKeyError, PyNotImplementedError, PyRuntimeError, PyTypeError, PyValueError,
};
use pyo3::prelude::*;
use pyo3::pyclass::IterNextOutput;
use pyo3::types::{PyIterator, PyList, PySequence, PyTuple};
use pyo3::FromPyObject;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::sync::{Arc, Mutex};
use std::time::Duration;
use std::{mem, thread};

/// Strand enum
#[pyclass]
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Strand {
    /// Maps to the forward strand
    Forward,
    /// Maps to the Reverse strand
    Reverse,
}

/// Enum containing results from multithreaded Alignment
#[derive(Debug, Clone)]
enum WorkQueue<T> {
    /// Lemme see you work work work work, shorty sumthin sumthin
    Work(T),
    /// The threads are finished
    Done,
    /// Result of multi threaded mapping queue
    Result(T),
    /// All threads have finished
    Finished,
}

/// Implement `Display` for `Strand`.
impl Display for Strand {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        let strand_string = match self {
            Strand::Forward => String::from("+"),
            Strand::Reverse => String::from("-"),
        };
        write!(f, "{strand_string}")
    }
}

#[pymethods]
impl Strand {
    /// String representation of the strand
    fn __str__(&self) -> String {
        format!("{}", &self)
    }
}

impl Strand {
    /// Convert the minimap2 crate strand to the strand enum found in our crate.
    fn from_mm2_strand(strand: minimap2::Strand) -> Strand {
        match strand {
            minimap2::Strand::Forward => Strand::Forward,
            minimap2::Strand::Reverse => Strand::Reverse,
        }
    }
}

/// Result of an alignment.
/// Attributes can be accessed using the `minimap2/mappy` attribute names, or
/// longer form methods.
///
/// # Examples
///
/// ```
///     use mappy_rs::{Mapping, Strand};
///     let m = Mapping {
///         query_start: 32,
///         query_end: 33,
///         strand: Strand::Forward,
///         target_name: String::from("Shark_bait"),
///         target_len: 10,
///         target_start: 10,
///         target_end: 11,
///         match_len: 10,
///         block_len: 10,
///         mapq: 69,
///         is_primary: true,
///         cigar: vec![(10, 11)],
///         NM: 10,
///         MD: None,
///         cs: None
///     };
///     // valid
///     assert!(m.target_start == 10); // also gets the mapping start
/// ```
///
/// In python this can be referenced as mapping.r_st or mapping.target_start.
#[pyclass]
#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(non_snake_case)]
pub struct Mapping {
    /// Mapping start on the query DNA sequence
    #[pyo3(get)]
    pub query_start: i32,
    /// Mapping end on the query DNA sequence
    #[pyo3(get)]
    pub query_end: i32,
    /// DNA strand the mapping is on (Forward or reverse)
    pub strand: Strand,
    /// The name of the chromosome/contig mapped to
    #[pyo3(get)]
    pub target_name: String,
    /// Length of the target contig mapped to
    #[pyo3(get)]
    pub target_len: i32,
    /// Mapping start on the contig mapped to
    #[pyo3(get)]
    pub target_start: i32,
    /// Mapping end on the contig mapped to
    #[pyo3(get)]
    pub target_end: i32,
    /// Match length of the alignment
    #[pyo3(get)]
    pub match_len: i32,
    /// Block length of the alignment (includes gaps)
    #[pyo3(get)]
    pub block_len: i32,
    /// Alignment quality between 0-255, with 255 for missing.
    #[pyo3(get)]
    pub mapq: u32,
    /// Alignment is primary or not
    #[pyo3(get)]
    pub is_primary: bool,
    /// The CIGAR operations/numbers of the alignment
    #[pyo3(get)]
    pub cigar: Vec<(u32, u8)>,
    /// Total number of matchs, mismatches and gaps in the alignment
    #[pyo3(get)]
    pub NM: i32,
    /// MD string of the alignment
    #[pyo3(get)]
    pub MD: Option<String>,
    /// CIGAR string
    #[pyo3(get)]
    pub cs: Option<String>,
}

/// Implement `Display` for `Mapping`. Writes out a paf formatted Mapping result.
/// NB. As the paf record spec describes, this will not include the `query name` and `query length`
/// fields.
impl Display for Mapping {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        let tp = if self.is_primary { "tp:A:P" } else { "tp:A:S" };
        let cigar = self.get_cigar_str().unwrap();
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
            self.query_start,
            self.query_end,
            self.strand,
            self.target_name,
            self.target_len,
            self.target_start,
            self.target_end,
            self.match_len,
            self.block_len,
            self.mapq,
            tp,
            cigar
        )
    }
}

/// Selection of python methods for an Alignment Mapping.
#[pymethods]
impl Mapping {
    /// Implement the python `repr()` method.
    fn __repr__(&self) -> String {
        format!("{self:#?}")
    }

    /// Implement the string representation in python.
    fn __str__(&self) -> String {
        format!("{self}")
    }

    /// Get the target name from a `Mapping`. Alias for `mappy.Alignment.ctg`
    #[getter(ctg)]
    fn get_target_name(&self) -> PyResult<String> {
        Ok(self.target_name.clone())
    }
    /// Get the target contig length from a `Mapping`. Alias for `mappy.Alignment.ctg_len`
    #[getter(ctg_len)]
    fn get_target_len(&self) -> PyResult<i32> {
        Ok(self.target_len)
    }

    /// Get the mapping start on the reference from a `Mapping`. Alias for `mappy.Alignment.r_st`
    #[getter(r_st)]
    fn get_target_start(&self) -> PyResult<i32> {
        Ok(self.target_start)
    }

    /// Get the mapping end on the reference from a `Mapping`. Alias for `mappy.Alignment.r_en`
    #[getter(r_en)]
    fn get_target_end(&self) -> PyResult<i32> {
        Ok(self.target_end)
    }

    /// Get the mapping start on the query sequence from a `Mapping`. Alias for `mappy.Alignment.q_st`
    #[getter(q_st)]
    fn get_query_start(&self) -> PyResult<i32> {
        Ok(self.query_start)
    }

    /// Get the mapping end on the query sequence from a `Mapping`. Alias for `mappy.Alignment.q_en`
    #[getter(q_en)]
    fn get_query_end(&self) -> PyResult<i32> {
        Ok(self.query_end)
    }

    /// Get the strand from a `Mapping`. Alias for `mappy.Alignment.strand`
    #[getter(strand)]
    fn get_strand(&self) -> PyResult<i32> {
        Ok(match self.strand {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        })
    }

    /// Get the alignment block length from a `Mapping`. Alias for `mappy.Alignment.blen`
    #[getter(blen)]
    fn get_block_len(&self) -> PyResult<i32> {
        Ok(self.block_len)
    }

    /// Get the match length from a `Mapping`. Alias for `mappy.Aignment.mlen`
    #[getter(mlen)]
    fn get_match_len(&self) -> PyResult<i32> {
        Ok(self.match_len)
    }

    /// Get the cigar string from a `Mapping`. Alias for `mappy.Alignment.cigar_str`
    #[getter(cigar_str)]
    fn get_cigar_str(&self) -> PyResult<String> {
        let strs = self
            .cigar
            .clone()
            .into_iter()
            .map(|(n, op)| {
                let c = match op {
                    0 => "M",
                    1 => "I",
                    2 => "D",
                    3 => "N",
                    4 => "S",
                    5 => "H",
                    6 => "P",
                    7 => "=",
                    8 => "X",
                    _ => return Err("Invalid CIGAR code `{op}`"),
                };
                Ok(format!("{n}{c}"))
            })
            .collect::<Result<Vec<_>, _>>();
        match strs {
            Ok(cstr) => Ok(cstr.join("")),
            Err(err) => Err(PyValueError::new_err(err)),
        }
    }

    /// Return whether this `Mapping` is a primary mapping. Alias for `mappy.Alignment.is_primary`
    #[getter(is_primary)]
    fn get_is_primary(&self) -> PyResult<bool> {
        Ok(self.is_primary)
    }
}

/// Aligner struct, mimicking minimap2's python interface
#[pyclass(unsendable)]
#[allow(clippy::type_complexity)]

pub struct Aligner {
    /// Inner minimap2::Aligner
    pub aligner: minimap2::Aligner,
    /// Number of mapping threads
    n_threads: usize,
    /// thread handles
    _handles: Arc<Mutex<Vec<std::thread::JoinHandle<()>>>>,
    /// stop the threads
    stop: Arc<Mutex<bool>>,
    /// Work queue stores strings to map and ids to get the corresponding dict back
    work_queue: Arc<ArrayQueue<WorkQueue<(usize, String)>>>,
    /// Results of the threads go here
    results_queue: Arc<ArrayQueue<WorkQueue<(Vec<Mapping>, usize)>>>,
}
// unsafe impl Send for Aligner {}

#[pymethods]
impl Aligner {
    /// Initialise a new Py Class Aligner
    /// Aligner struct, mimicking minimap2's python interface
    #[new]
    #[pyo3(signature = (fn_idx_in=None, preset=None, k=None, w=None, min_cnt=None, min_chain_score=None, min_dp_score=None, bw=None, best_n=None, n_threads=3, fn_idx_out=None, max_frag_len=None, extra_flags=None, seq=None, scoring=None))]
    #[allow(clippy::too_many_arguments, unused_assignments)]
    fn py_new(
        fn_idx_in: Option<std::path::PathBuf>,
        preset: Option<String>,
        k: Option<usize>,
        w: Option<usize>,
        min_cnt: Option<usize>,
        min_chain_score: Option<usize>,
        min_dp_score: Option<usize>,
        bw: Option<usize>,
        best_n: Option<usize>,
        n_threads: usize,
        fn_idx_out: Option<std::path::PathBuf>,
        max_frag_len: Option<usize>,
        extra_flags: Option<usize>,
        seq: Option<String>,
        scoring: Option<&PyTuple>,
    ) -> PyResult<Self> {
        let mut mapopts = minimap2::MapOpt::default();
        let mut idxopts = minimap2::IdxOpt::default();
        unsafe { minimap2_sys::mm_set_opt(std::ptr::null(), &mut idxopts, &mut mapopts) };
        if let Some(preset) = preset {
            let _preset = std::ffi::CString::new(preset).unwrap();
            unsafe { minimap2_sys::mm_set_opt(_preset.as_ptr(), &mut idxopts, &mut mapopts) };
        }
        // For 'drop-in' mappy compatibility we should add the flag 4
        mapopts.flag |= 4;
        idxopts.batch_size |= 0x7fffffffffffffff_u64;

        if let Some(k) = k {
            idxopts.k = k as i16
        }
        if let Some(w) = w {
            idxopts.w = w as i16
        }
        if let Some(min_cnt) = min_cnt {
            mapopts.min_cnt = min_cnt as i32
        }
        if let Some(min_chain_score) = min_chain_score {
            mapopts.min_chain_score = min_chain_score as i32
        }
        if let Some(min_dp_score) = min_dp_score {
            mapopts.min_dp_max = min_dp_score as i32
        }
        if let Some(bw) = bw {
            mapopts.bw = bw as i32
        }
        if let Some(best_n) = best_n {
            mapopts.best_n = best_n as i32
        }
        if let Some(max_frag_len) = max_frag_len {
            mapopts.max_frag_len = max_frag_len as i32
        }
        if let Some(extra_flags) = extra_flags {
            mapopts.flag |= extra_flags as i64
        }
        if let Some(scoring) = scoring {
            if scoring.len() >= 4 {
                mapopts.a = scoring.get_item(0).unwrap().extract::<i32>().unwrap();
                mapopts.b = scoring.get_item(1).unwrap().extract::<i32>().unwrap();
                mapopts.q = scoring.get_item(2).unwrap().extract::<i32>().unwrap();
                mapopts.e = scoring.get_item(3).unwrap().extract::<i32>().unwrap();
                mapopts.q2 = mapopts.q;
                mapopts.e2 = mapopts.e;
                if scoring.len() >= 6 {
                    mapopts.q2 = scoring.get_item(4).unwrap().extract::<i32>().unwrap();
                    mapopts.e2 = scoring.get_item(5).unwrap().extract::<i32>().unwrap();
                    if scoring.len() >= 7 {
                        mapopts.sc_ambi = scoring.get_item(6).unwrap().extract::<i32>().unwrap();
                    }
                }
            }
        }

        // TODO: The scoping rules are tricky here - maybe
        if let Some(_seq) = seq {
            return Err(PyNotImplementedError::new_err("Not Implemented"));
        }
        if let Some(_fn_idx_out) = fn_idx_out {
            // If this is set, we create an MMI but cannot use it
            return Err(PyNotImplementedError::new_err("Not Implemented"));
        }
        if let Some(fn_idx_in) = fn_idx_in {
            let fn_in = std::ffi::CString::new(fn_idx_in.to_str().unwrap()).unwrap();
            let idx_reader = std::mem::MaybeUninit::new(unsafe {
                minimap2_sys::mm_idx_reader_open(fn_in.as_ptr(), &idxopts, std::ptr::null())
            });

            let mut idx: std::mem::MaybeUninit<*mut minimap2_sys::mm_idx_t> =
                std::mem::MaybeUninit::uninit();

            let idx_reader = unsafe { idx_reader.assume_init() };

            unsafe {
                idx = std::mem::MaybeUninit::new(minimap2_sys::mm_idx_reader_read(
                    &mut *idx_reader as *mut minimap2_sys::mm_idx_reader_t,
                    n_threads as libc::c_int,
                ));
                // Close the reader
                minimap2_sys::mm_idx_reader_close(idx_reader);
                // Set index opts
                minimap2_sys::mm_mapopt_update(&mut mapopts, *idx.as_ptr());
                // Idx index name
                minimap2_sys::mm_idx_index_name(idx.assume_init());
            };
            let al = Aligner {
                aligner: minimap2::Aligner {
                    mapopt: mapopts,
                    idxopt: idxopts,
                    threads: n_threads,
                    idx: Some(unsafe { *idx.assume_init() }),
                    idx_reader: Some(unsafe { *idx_reader }),
                },
                n_threads: 0,
                _handles: Arc::new(Mutex::new(vec![])),
                stop: Arc::new(Mutex::new(false)),
                work_queue: Arc::new(ArrayQueue::<WorkQueue<(usize, String)>>::new(50000)),
                results_queue: Arc::new(ArrayQueue::<WorkQueue<(Vec<Mapping>, usize)>>::new(50000)),
            };
            // al.setup_signal();
            return Ok(al);
        }
        Err(PyRuntimeError::new_err("Did not create or open an index"))
    }

    /// Return the sequence names contained within an index as a list.
    #[getter]
    fn seq_names(&self) -> PyResult<Vec<String>> {
        if !self.aligner.has_index() {
            return Err(PyRuntimeError::new_err("Index hasn't loaded"));
        }
        unsafe {
            let ns = (self.aligner.idx.unwrap()).n_seq;
            let mut sn = Vec::with_capacity(ns as usize);
            for i in 0..ns {
                sn.push(
                    std::ffi::CStr::from_ptr(
                        (*((self.aligner.idx.unwrap()).seq.offset(i as isize))).name,
                    )
                    .to_str()
                    .unwrap()
                    .to_string(),
                )
            }
            Ok(sn)
        }
    }

    ///  Retrieves a (sub)sequence from the index and returns it as a Python string. None is
    ///  returned if name is not present in the index or the start/end coordinates are invalid
    ///  or if the index does not contain any sequence.
    #[pyo3(signature = (name, start=0, end=2147483647), text_signature = "(name, start=0, end=2147483647)")]
    fn seq(&self, name: String, start: i32, end: i32) -> PyResult<Option<String>> {
        Ok(match self._get_index_seq(name, start, end) {
            Ok(res) => Some(res),
            Err(_) => None,
        })
    }

    /// Map a single read, blocking
    #[pyo3(signature = (seq, seq2=None, cs=false, MD=false), text_signature = "(seq, seq2=None, cs=False, MD=False)")]
    #[allow(non_snake_case)]
    fn map(&self, seq: String, seq2: Option<String>, cs: bool, MD: bool) -> PyResult<Vec<Mapping>> {
        // TODO: PyIterProtocol to map single reads and return as a generator
        if let Some(_seq2) = seq2 {
            return Err(PyNotImplementedError::new_err(
                "Using `seq2` is not implemented",
            ));
        }
        match self.aligner.map(
            seq.as_bytes(),
            cs,
            MD,
            Some(self.aligner.mapopt.max_frag_len as usize),
            None,
        ) {
            Ok(res) => Ok(res
                .into_iter()
                .map(|m| {
                    let a = m.alignment.unwrap();
                    Mapping {
                        query_start: m.query_start,                // i32,
                        query_end: m.query_end,                    // i32,
                        strand: Strand::from_mm2_strand(m.strand), // Strand,
                        target_name: m.target_name.unwrap(),       // String,
                        target_len: m.target_len,                  // i32,
                        target_start: m.target_start,              // i32,
                        target_end: m.target_end,                  // i32,
                        match_len: m.match_len,                    // i32,
                        block_len: m.block_len,                    // i32,
                        mapq: m.mapq,                              // u32,
                        is_primary: m.is_primary,                  // bool
                        cigar: a.cigar.unwrap_or_default(),        // Vec<(u32, u8)>
                        NM: a.nm,                                  // i32
                        MD: a.md,                                  // Option<String>
                        cs: a.cs,                                  // Option<String>
                    }
                })
                .collect()),
            Err(e) => Err(PyRuntimeError::new_err(e)),
        }
    }

    /// Map a single read, blocking
    #[pyo3(signature = (_seq, seq2=None, _cs=false, _MD=false), text_signature = "(_seq, seq2=None, _cs=False, _MD=False)")]
    #[allow(non_snake_case)]
    fn map_no_op(
        &self,
        _seq: String,
        seq2: Option<String>,
        _cs: bool,
        _MD: bool,
    ) -> PyResult<Vec<Mapping>> {
        // TODO: PyIterProtocol to map single reads and return as a generator
        if let Some(_seq2) = seq2 {
            return Err(PyNotImplementedError::new_err(
                "Using `seq2` is not implemented",
            ));
        }
        Ok(self.no_op_map())
    }

    ///  Enable multi threading on this mappy instance.
    ///
    /// Example
    /// -------
    /// `aligner::enable_threading(8)`
    #[pyo3(signature = (n_threads), text_signature = "(n_threads=8)")]
    fn enable_threading(&mut self, n_threads: usize) -> PyResult<()> {
        self.n_threads = n_threads;
        let dones = Arc::new(Mutex::new(vec![false; n_threads]));
        for i in 0..n_threads {
            let _aligner = self.aligner.clone();
            let stop = Arc::clone(&self.stop);
            let wq = Arc::clone(&self.work_queue);
            let rq = Arc::clone(&self.results_queue);
            let thread_number = i;
            let done_ref = Arc::clone(&dones);

            // start the threads
            std::thread::spawn(move || {
                loop {
                    // STOP SIGNAL RECEVIED SIGINT/SIGTERM
                    if *stop.lock().unwrap() {
                        break;
                    }
                    let wait_as_done = {
                        let mut dr: std::sync::MutexGuard<'_, Vec<bool>> = done_ref.lock().unwrap();
                        let done = *dr.get(thread_number).unwrap();
                        let all_done = all(dr.iter(), |elt| *elt);
                        if all_done {
                            // everythread has sent a done, so set them all to not done again
                            for b in dr.iter_mut() {
                                *b = !*b;
                            }
                        }
                        done & !all_done
                    };
                    // this thread has sent a done and not all other threads are fininshed
                    if wait_as_done {
                        std::thread::sleep(Duration::from_millis(1));
                        continue;
                    }
                    match wq.pop() {
                        None => std::thread::sleep(Duration::from_millis(10)),
                        Some(work_item) => {
                            match work_item {
                                WorkQueue::Done => {
                                    rq.push(WorkQueue::Done).unwrap();
                                    {
                                        done_ref.lock().unwrap()[thread_number] = true;
                                    }
                                }
                                WorkQueue::Work((id_num, seq)) => {
                                    match _aligner.map(
                                        seq.as_bytes(),
                                        true,
                                        false,
                                        Some(_aligner.mapopt.max_frag_len as usize),
                                        None,
                                    ) {
                                        Ok(_mappings) => {
                                            mem::drop(seq);
                                            let mappings: Vec<Mapping> = _mappings
                                                .into_iter()
                                                .map(|m| {
                                                    let a = m.alignment.unwrap();
                                                    Mapping {
                                                        query_start: m.query_start,                // i32,
                                                        query_end: m.query_end, // i32,
                                                        strand: Strand::from_mm2_strand(m.strand), // Strand,
                                                        target_name: m.target_name.unwrap(), // String,
                                                        target_len: m.target_len,            // i32,
                                                        target_start: m.target_start,        // i32,
                                                        target_end: m.target_end,            // i32,
                                                        match_len: m.match_len,              // i32,
                                                        block_len: m.block_len,              // i32,
                                                        mapq: m.mapq,                        // u32,
                                                        is_primary: m.is_primary,            // bool
                                                        cigar: a.cigar.unwrap_or_default(), // Vec<(u32, u8)>
                                                        NM: a.nm,
                                                        MD: a.md,
                                                        cs: a.cs,
                                                    }
                                                })
                                                .collect();
                                            rq.push(WorkQueue::Result((mappings, id_num))).unwrap();
                                        }
                                        Err(_) => {
                                            eprintln!("Failed to map sequence in threaded implementation.")
                                        }
                                    }
                                }
                                _ => {
                                    println!("What is this doing in the work queue")
                                }
                            }
                        }
                    }
                }
            });
        }
        Ok(())
    }

    /// Align a sequence Optionally back off if we fail to add the sequence to the queue, in the case that the work queue is full.
    #[pyo3(signature = (seqs, back_off=true))]
    fn map_batch(&self, seqs: &PyAny, back_off: bool) -> PyResult<AlignmentBatchResultIter> {
        let mut res = AlignmentBatchResultIter::new();
        // Set the number of threads
        res.set_n_threads(self.n_threads);
        // do the heavy work
        self._map_batch(&mut res, seqs, back_off)?;
        // let return_metadata: (i32, i32, String) = (metadata.read_number, metadata.channel_number, String::from("hdea"));
        Ok(res)
    }

    /// Return whether or not this Aligner has an index.
    fn __bool__(&self) -> PyResult<bool> {
        Ok(self.aligner.idx.is_some())
    }

    /// Get the k value from the index.
    #[getter]
    fn k(&self) -> PyResult<i32> {
        Ok(self.aligner.idx.unwrap().k)
    }
    /// Get the w value form the index.
    #[getter]
    fn w(&self) -> PyResult<i32> {
        Ok(self.aligner.idx.unwrap().w)
    }

    /// Get the number of sequences present in the index
    #[getter]
    fn n_seq(&self) -> PyResult<u32> {
        Ok(self.aligner.idx.unwrap().n_seq)
    }
}

impl Aligner {
    /// Instead of calling out to the ALigner, return a predefined dummy mapping
    pub fn no_op_map(&self) -> Vec<Mapping> {
        vec![Mapping {
            query_start: 0,                                             // i32,
            query_end: 1000,                                            // i32,
            strand: Strand::from_mm2_strand(minimap2::Strand::Forward), // Strand,
            target_name: String::from("Hello"),                         // String,
            target_len: 101010,                                         // i32,
            target_start: 10,                                           // i32,
            target_end: 1010,                                           // i32,
            match_len: 1000,                                            // i32,
            block_len: 1000,                                            // i32,
            mapq: 60,                                                   // u32,
            is_primary: true,                                           // bool
            cigar: vec![],                                              // Vec<(u32, u8)>
            NM: 0,
            MD: None,
            cs: Some(String::from("Cigar string")),
        }]
    }
    /// Setup signal catching for ctrl c to stop threads
    pub fn setup_signal(&self) {
        let stop = Arc::clone(&self.stop);
        ctrlc::set_handler(move || {
            println!("Signal intercepted");
            *stop.lock().unwrap() = true;
            std::process::exit(0);
        })
        .expect("Failed to set signal listener");
    }
    /// Private function
    /// Get a sequence or subsequence of a contig loaded into the index.
    pub fn _get_index_seq(&self, name: String, start: i32, mut end: i32) -> Result<String, &str> {
        if !self.aligner.has_index() {
            return Err("No index");
        }
        if (self.aligner.mapopt.flag & minimap2_sys::MM_F_CIGAR as i64 != 0)
            && (self.aligner.idx.unwrap().flag & minimap2_sys::MM_I_NO_SEQ as i32 != 0)
        {
            return Err("No sequence in this index");
        }
        let ref_seq_id: i32 = unsafe {
            minimap2_sys::mm_idx_name2id(
                self.aligner.idx.as_ref().unwrap() as *const minimap2_sys::mm_idx_t,
                std::ffi::CString::new(name)
                    .unwrap()
                    .as_bytes_with_nul()
                    .as_ptr() as *const i8,
            )
        };
        if (ref_seq_id < 0) | (ref_seq_id as u32 >= self.aligner.idx.unwrap().n_seq) {
            return Err("Could not find reference in index");
        }

        let ref_seq_offset = unsafe {
            *(self
                .aligner
                .idx
                .as_ref()
                .unwrap()
                .seq
                .offset(ref_seq_id as isize))
        };
        let ref_seq_len = ref_seq_offset.len as i32;
        if start >= ref_seq_len || start >= end {
            return Err("Funky start and end coords");
        }
        if end < 0 || end > ref_seq_len {
            end = ref_seq_len;
        }
        let seq_len = end - start;
        let mut seq_buf: Vec<u8> = vec![0; seq_len as usize];
        let _len = unsafe {
            minimap2_sys::mm_idx_getseq(
                self.aligner.idx.as_ref().unwrap() as *const minimap2_sys::mm_idx_t,
                ref_seq_id as u32,
                start as u32,
                end as u32,
                seq_buf.as_mut_ptr(),
            )
        };
        for c in &mut seq_buf {
            *c = match *c {
                0 => 65, // A
                1 => 67, // C
                2 => 71, // G
                3 => 84, // T
                4 => 78, // N
                _ => return Err("Got an unknown char, not {ACGTN}"),
            }
        }
        Ok(std::string::String::from_utf8(seq_buf).unwrap())
    }

    /// Align a batch of reads provided in an iterator, using a threadpool with the number of threads specified by
    /// .enable_threading()
    #[allow(clippy::type_complexity)]
    pub fn _map_batch(
        &self,
        res: &mut AlignmentBatchResultIter,
        seqs: &PyAny,
        back_off: bool,
    ) -> PyResult<()> {
        if self.n_threads == 0_usize {
            return Err(PyRuntimeError::new_err(
                "Multi threading not enabled on this instance. Please call `.enable_threading()`",
            ));
        }
        match seqs.extract() {
            Ok(SupportedTypes::List(_)) => (),
            Ok(SupportedTypes::Tuple(_)) => (),
            Ok(SupportedTypes::Iter(_)) => (),
            Ok(SupportedTypes::Sequence(_)) => (),
            _ => {
                return Err(PyTypeError::new_err(
                    "Unsupported batch type, pass a list, iter, generator or tuple",
                ))
            }
        };
        let results_queue: Arc<ArrayQueue<WorkQueue<(Vec<Mapping>, usize)>>> =
            Arc::clone(&self.results_queue);
        let results_tx = res.tx.clone();
        let counter = Arc::clone(&res._n_finished_threads);
        let n_threads = res._n_threads;
        std::thread::spawn(move || {
            loop {
                //             // pop returns None if the queue is empty, which is possible at the start as data hasn't been added below
                match results_queue.pop() {
                    //                 // We
                    Some(worky) => match worky {
                        // each thread can only see one workqueue DONE
                        WorkQueue::Done => {
                            // This thread has finished
                            // Lock the mutex and increment
                            let mut num = counter.lock().unwrap();
                            *num += 1;
                            // println!("{num}");
                            // ALL threads have finished
                            if *num == n_threads {
                                results_tx.send(WorkQueue::Finished).unwrap();
                                // reset number of finshed threads
                                break;
                            }
                        }
                        WorkQueue::Result(result) => {
                            let id = result.1;
                            match results_tx.send(WorkQueue::Result(result)) {
                                Ok(()) => {}
                                Err(e) => {
                                    eprintln!("Internal error returning data, the receiver iterator has finished. {e} {id}");
                                    break;
                                }
                            }
                        }
                        _ => {
                            eprintln!("Wrong WorkQueue arm seen in worker thread.")
                        }
                    },
                    // Todo Crossbeam backoff rather than continue
                    None => {
                        std::thread::sleep(Duration::from_millis(5));
                    }
                }
                //             // (id_num, seq): (usize, String)
            }
        });
        let iter = match seqs.iter() {
            Ok(it) => it,
            _ => return Err(PyTypeError::new_err("Could not iterate batch")),
        };
        let work_queue: Arc<ArrayQueue<WorkQueue<(usize, String)>>> = Arc::clone(&self.work_queue);
        for (id_num, py_dict) in iter.enumerate() {
            let py_dict = py_dict?;
            let data: HashMap<String, Py<PyAny>> = match py_dict.extract() {
                Ok(x) => x,
                _ => {
                    return Err(PyTypeError::new_err(
                        "Element in iterable is not a dictionary",
                    ))
                }
            };
            res.data.insert(id_num, data);
            let seq: String = match py_dict.get_item("seq") {
                Ok(seq) => match seq.extract::<String>() {
                    Ok(seq) => seq.clone(),
                    _ => return Err(PyValueError::new_err("`seq` must be a string")),
                },
                _ => {
                    return Err(PyKeyError::new_err(
                        "AHHH Key 🗝️  not found in iterated dictionary",
                    ))
                }
            };
            match work_queue.push(WorkQueue::Work((id_num, seq))) {
                Ok(()) => {}
                Err(e) => {
                    if back_off {
                        let mut attempts = 0;
                        let max_attempts = 6;
                        let mut sleep_duration = Duration::from_millis(50); // Initial sleep duration (in milliseconds)

                        while attempts < max_attempts {
                            if work_queue.push(e.clone()).is_ok() {
                                break; // Operation succeeded
                            }

                            attempts += 1;
                            thread::sleep(sleep_duration);

                            // Increase the sleep duration exponentially
                            sleep_duration *= 2;
                        }
                        if attempts == 6 {
                            eprintln!("Internal error adding data to work queue, with backoff. {e:#?}, {id_num}, Attempts: {attempts}");
                        }
                    } else {
                        eprintln!("Internal error adding data to work queue, without backoff. {e:#?} {id_num}");
                        return Err(PyErr::new::<PyRuntimeError, _>(format!(
                            "Internal error adding data to work queue, without backoff. {e:#?} {id_num}. Is your fastq batch larger than 50000? Perhaps try `map_batch` with back_off=True?",
                            e = e,
                            id_num = id_num
                        )));
                    }
                }
            }
        }
        // Now we add n_thread dones, one for each thread. When the threads see this they know to close as there is no more data
        for _ in 0..self.n_threads {
            work_queue.push(WorkQueue::Done).unwrap();
        }

        Ok(())
    }
}

/// Python iterable types that are accepted by the `Aligner.map_batch()` function
#[derive(FromPyObject)]
enum SupportedTypes<'py> {
    /// PyList wrapper
    List(&'py PyList),
    /// PyTuple wrapper
    Tuple(&'py PyTuple),
    /// PyIter wrapper
    Iter(&'py PyIterator),
    /// PySequence wrapper
    Sequence(&'py PySequence),
}

/// Struct for returning data to the python runtime as an iterabled.
#[pyclass]
pub struct AlignmentBatchResultIter {
    /// Sender of results into this scope
    tx: Sender<WorkQueue<(Vec<Mapping>, usize)>>,
    /// Receive the sent data
    rx: Receiver<WorkQueue<(Vec<Mapping>, usize)>>,
    /// HashMap for caching sent data
    data: FnvHashMap<usize, HashMap<String, Py<PyAny>>>,
    /// Number of threads, which checks against the number offinished threads
    _n_threads: usize,
    /// Number of finished threads, used to know when to close the receiver. Is unlocked in the worker threads.
    _n_finished_threads: Arc<Mutex<usize>>,
}

impl Default for AlignmentBatchResultIter {
    /// Impl default for the `AlignmentBatchResultIter.`
    fn default() -> Self {
        Self::new()
    }
}

/// Iterator for the batch results from a multi threaded call to mapper
#[pymethods]
impl AlignmentBatchResultIter {
    /// Initialise a new `AlignmentBatchResultIter`. Spawns the Send And Receive channels.
    #[new]
    pub fn new() -> Self {
        let (tx, rx) = bounded(20000);
        AlignmentBatchResultIter {
            tx,
            rx,
            data: FnvHashMap::default(),
            _n_threads: 0_usize,
            _n_finished_threads: Arc::new(Mutex::new(0_usize)),
        }
    }

    /// Set the number of threads on this mappy_rs instance
    pub fn set_n_threads(&mut self, n_threads: usize) {
        self._n_threads = n_threads;
    }

    /// Returns the Iterable, in this case the struct itself.
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    /// Returns the next element in the Iterator.
    #[allow(clippy::type_complexity)]
    fn __next__(&mut self) -> IterNextOutput<(Vec<Mapping>, HashMap<String, Py<PyAny>>), &str> {
        let try_recv = self.rx.recv();
        match try_recv {
            Ok(work_queue_member) => match work_queue_member {
                WorkQueue::Finished => IterNextOutput::Return("Finished"),
                WorkQueue::Result((mapping, id_num)) => {
                    let data = self.data.remove(&id_num).unwrap();
                    IterNextOutput::Yield((mapping, data))
                }
                _ => {
                    eprintln!("Received wrong variant as a Result");
                    IterNextOutput::Return("Wrong variant")
                }
            },
            Err(RecvError) => {
                eprintln!("Receiver Error");
                IterNextOutput::Return("Receiver error - channel was closed")
            }
        }
    }
}

/// Initialise the python module and add the Aligner class.
#[pymodule]
fn mappy_rs(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<Aligner>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn get_resource_dir() -> PathBuf {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("resources/test");
        path
    }

    fn get_test_file(file: &str) -> PathBuf {
        let mut path = get_resource_dir();
        path.push(file);
        path
    }

    fn get_test_aligner() -> Result<Aligner, PyErr> {
        let path = get_test_file("test.mmi");
        Aligner::py_new(
            Some(path),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            4_usize,
            None,
            None,
            None,
            None,
            None,
        )
    }

    #[test]
    fn load_index() {
        let al = get_test_aligner().unwrap();
        assert!(al.aligner.has_index());
    }

    #[test]
    fn test_property_k() {
        let al = get_test_aligner().unwrap();
        assert!(al.k().unwrap() == 15);
    }

    #[test]
    fn test_property_w() {
        let al = get_test_aligner().unwrap();
        assert!(al.w().unwrap() == 10);
    }

    #[test]
    fn test_property_n_seq() {
        let al = get_test_aligner().unwrap();
        assert!(al.n_seq().unwrap() == 4);
    }

    #[test]
    fn test_property_seq_names() {
        let al = get_test_aligner().unwrap();
        let expected = vec![
            "Bacillus_subtilis",
            "Enterococcus_faecalis",
            "Escherichia_coli_1",
            "Escherichia_coli_2",
        ];
        let mut seq_names = al.seq_names().unwrap();
        seq_names.sort();
        assert!(seq_names == expected);
    }

    #[test]
    fn test_get_seq() {
        let al = get_test_aligner().unwrap();
        let contig = "Bacillus_subtilis";
        let expected = "AGAGTGAAGCCAATATTCCGATAACGATTGCTTTCATGATATCCCTCATTCTGGCATTATTTTTTTATACTATACTATTC\
                        GATATCGCACAGATCAATGGAGTCGTGAGAAAATAAACATGTTTTGCGAACCGCTATGTGTGGAAGACAAAAAATGGAGG\
                        TGAAATTGATGGAAGCAAAGACACAGGCGTACTTTTTTCAGGATGATGGCAGGATTCCGAATCACCCTGATTTTCCGCTC\
                        GTTGTGTATCAAAACGCACTCAAGGACACCGGTCAGGCAGAGCGGATCGTCAACCGGCATGGCTGGTCAAACAGCTGGTC\
                        GGGGAGTGTTTTTCCATACCATCATTATCACAGCAATACGCATGAAGTCCTGATTGCAGTTCGGGGAGAGGCTGTGATTC";
        let seq = al
            .seq(String::from(contig), 0, 2147483647)
            .unwrap()
            .unwrap();
        assert!(seq == *expected);
    }

    #[test]
    fn map_one() {
        let al = get_test_aligner().unwrap();
        let mappings = al.map(
            String::from("AGAGCAGGTAGGATCGTTGAAAAAAGAGTACTCAGGATTCCATTCAACTTTTACTGATTTGAAGCGTACTGTTTATGGCC\
                          AAGAATATTTACGTCTTTACAACCAATACGCAAAAAAAGGTTCATTGAGTTTGGTTGTGATTTGATGAAAATTACTGAGA\
                          ATAACAGGATTATTAAGCTGATTGATGAACTAAATCAGCTTAATAAATATTCTTTGCAGATAGGAATATTTGGGGAAAAT\
                          GATTCTTTTATGGCGATGTTGGCCCAAGTTCATGAATTTGGGGTGACTATTCGTCCCAAAGGTCGTTTTCTTGTTATACC\
                          ACTTATGAAAAAGTATAGAGGTAAAAGTCCACGTCAATTTGATTTGTTTTTTATGCAAACTAAAGAAAATCACAAGTTTT"),
            None, true, false).unwrap();
        assert!(mappings.len() == 1);
        assert!(mappings[0].get_target_start().unwrap() == 0);
        assert!(mappings[0].get_target_end().unwrap() == 400);
    }
}
