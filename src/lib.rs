use std::cell::RefCell;
use std::mem;
use std::mem::MaybeUninit;
use std::path::Path;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::Mutex;
use std::time::Duration;

use crossbeam::queue::ArrayQueue;
use minimap2_sys::*;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::pyclass::IterNextOutput;

// Pyiterator for handling generators passed in
use pyo3::types::{PyIterator, PyTuple};
use simdutf8::basic::from_utf8;

/// Alias for mm_mapop_t
pub type MapOpt = mm_mapopt_t;

/// Alias for mm_idxopt_t
pub type IdxOpt = mm_idxopt_t;

// TODO: Probably a better way to handle this...
static MAP_ONT: &str = "map-ont\0";
static AVA_ONT: &str = "ava-ont\0";
static MAP10K: &str = "map10k\0";
static AVA_PB: &str = "ava-pb\0";
static MAP_HIFI: &str = "map-hifi\0";
static ASM: &str = "asm\0";
static ASM5: &str = "asm5\0";
static ASM10: &str = "asm10\0";
static ASM20: &str = "asm20\0";
static SHORT: &str = "short\0";
static SR: &str = "sr\0";
static SPLICE: &str = "splice\0";
static CDNA: &str = "cdna\0";

#[derive(Debug)]
enum WorkQueue<T> {
    Work(T),
    Done,
    Starting,
}
/// Strand enum
#[pyclass]
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Strand {
    Forward,
    Reverse,
}

/// Preset's for minimap2 config
#[derive(Debug)]
pub enum Preset {
    MapOnt,
    AvaOnt,
    Map10k,
    AvaPb,
    MapHifi,
    Asm,
    Asm5,
    Asm10,
    Asm20,
    Short,
    Sr,
    Splice,
    Cdna,
}

// Convert to c string for input into minimap2
impl From<Preset> for *const i8 {
    fn from(preset: Preset) -> Self {
        match preset {
            Preset::MapOnt => MAP_ONT.as_bytes().as_ptr() as *const i8,
            Preset::AvaOnt => AVA_ONT.as_bytes().as_ptr() as *const i8,
            Preset::Map10k => MAP10K.as_bytes().as_ptr() as *const i8,
            Preset::AvaPb => AVA_PB.as_bytes().as_ptr() as *const i8,
            Preset::MapHifi => MAP_HIFI.as_bytes().as_ptr() as *const i8,
            Preset::Asm => ASM.as_bytes().as_ptr() as *const i8,
            Preset::Asm5 => ASM5.as_bytes().as_ptr() as *const i8,
            Preset::Asm10 => ASM10.as_bytes().as_ptr() as *const i8,
            Preset::Asm20 => ASM20.as_bytes().as_ptr() as *const i8,
            Preset::Short => SHORT.as_bytes().as_ptr() as *const i8,
            Preset::Sr => SR.as_bytes().as_ptr() as *const i8,
            Preset::Splice => SPLICE.as_bytes().as_ptr() as *const i8,
            Preset::Cdna => CDNA.as_bytes().as_ptr() as *const i8,
        }
    }
}

/// Alignment type
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlignmentType {
    Primary,
    Secondary,
    Inversion,
}

/// Alignment struct when alignment flag is set
#[pyclass]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub is_primary: bool,
    pub cigar: Option<String>,
}

#[pymethods]
impl Alignment {
    // For `__repr__` we want to return a string that Python code could use to recreate
    // the `Alignment`.
    fn __repr__(&self) -> String {
        // We use the `format!` macro to create a string. Its first argument is a
        // format string, followed by any number of parameters which replace the
        // `{}`'s in the format string.
        //
        //                       👇 Tuple field access in Rust uses a dot
        format!("{:#?}", self)
    }
}

/// Mapping result
#[pyclass]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mapping {
    // The query sequence name.
    #[pyo3(get)]
    pub query_name: Option<String>,
    #[pyo3(get)]
    pub query_len: Option<u32>,
    #[pyo3(get)]
    pub query_start: i32,
    #[pyo3(get)]
    pub query_end: i32,
    #[pyo3(get)]
    pub strand: Strand,
    #[pyo3(get)]
    pub target_name: Option<String>,
    #[pyo3(get)]
    pub target_len: i32,
    #[pyo3(get)]
    pub target_start: i32,
    #[pyo3(get)]
    pub target_end: i32,
    #[pyo3(get)]
    pub match_len: i32,
    #[pyo3(get)]
    pub block_len: i32,
    #[pyo3(get)]
    pub mapq: u32,
    #[pyo3(get)]
    pub alignment: Option<Alignment>,
}

#[pymethods]
impl Mapping {
    // #[getter]
    // fn name(&self) -> PyResult<String> {
    //     Ok(self.query_name.as_ref().unwrap().to_string())
    // }

    // #[setter]
    // fn set_name(&mut self, value: String) -> PyResult<()> {
    //     self.query_name = Some(value);
    //     Ok(())
    // }

    // For `__repr__` we want to return a string that Python code could use to recreate
    // the `Number`, like `Number(5)` for example.
    fn __repr__(&self) -> String {
        // We use the `format!` macro to create a string. Its first argument is a
        // format string, followed by any number of parameters which replace the
        // `{}`'s in the format string.
        //
        //                       👇 Tuple field access in Rust uses a dot
        format!("{:#?}", self)
    }

    // `__str__` is generally used to create an "informal" representation, so we
    // just forward to `i32`'s `ToString` trait implementation to print a bare number.
    // fn __str__(&self) -> String {
    //     "Get mapped scrub".to_string()
    // }
    // todo return paf formatted results
}

// Thread local buffer (memory management) for minimap2
thread_local! {
    static BUF: RefCell<ThreadLocalBuffer> = RefCell::new(ThreadLocalBuffer::new());
}

/// ThreadLocalBuffer for minimap2 memory management
struct ThreadLocalBuffer {
    buf: *mut mm_tbuf_t,
}

impl ThreadLocalBuffer {
    pub fn new() -> Self {
        let buf = unsafe { mm_tbuf_init() };
        Self { buf }
    }
}

/// Handle destruction of thread local buffer properly.
impl Drop for ThreadLocalBuffer {
    fn drop(&mut self) {
        unsafe { mm_tbuf_destroy(self.buf) };
    }
}

impl Default for ThreadLocalBuffer {
    fn default() -> Self {
        Self::new()
    }
}

#[pyclass]
#[derive(Debug)]
struct AlignmentResult {
    mappings: std::vec::IntoIter<Mapping>,
    #[pyo3(get)]
    metadata: (i32, i32, String),
}

#[pymethods]
impl AlignmentResult {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Mapping> {
        slf.mappings.next()
    }
}
/// Aligner struct, mimicking minimap2's python interface
///
/// ```
/// # use mappy_rs::*;
/// Aligner::builder();
/// ```
#[pyclass]
#[derive(Clone)]
pub struct Aligner {
    /// Index options passed to minimap2 (mm_idxopt_t)
    pub idxopt: IdxOpt,

    /// Mapping options passed to minimap2 (mm_mapopt_t)
    pub mapopt: MapOpt,

    /// Number of threads to create the index with
    pub idx_threads: usize,

    // number of threads to align with
    pub map_threads: usize,

    /// Index created by minimap2
    pub idx: Option<mm_idx_t>,

    /// Index reader created by minimap2
    pub idx_reader: Option<mm_idx_reader_t>,

    work_queue: Arc<ArrayQueue<WorkQueue<String>>>,
}
unsafe impl Send for Aligner {}

/// Create a default aligner
impl Default for Aligner {
    fn default() -> Self {
        Self {
            idxopt: Default::default(),
            mapopt: Default::default(),
            idx_threads: 1,
            map_threads: 1,
            idx: None,
            idx_reader: None,
            work_queue: Arc::new(ArrayQueue::<WorkQueue<String>>::new(50000)),
        }
    }
}

impl Aligner {
    /// Create a new server builder that can configure a [`Server`].
    pub fn builder() -> Self {
        Aligner {
            mapopt: MapOpt {
                ..Default::default()
            },
            ..Default::default()
        }
    }
}

struct MetaData {
    read_id: String,
    read_number: i32,
    channel_number: i32,
}

impl From<(i32, i32, String)> for MetaData {
    fn from(item: (i32, i32, String)) -> Self {
        MetaData {
            read_id: item.2,
            read_number: item.0,
            channel_number: item.1,
        }
    }
}

// Implement python methods
// yay
#[pymethods]
impl Aligner {
    #[new]
    fn py_new(n_threads: usize, index: PathBuf) -> PyResult<Self> {
        if n_threads == 0 {
            Err(PyValueError::new_err("n_threads cannot be zero"))
        } else {
            let mut aligner = Aligner::builder()
                .preset(Preset::MapOnt)
                .with_idx_threads(2)
                .with_map_threads(n_threads)
                .with_cigar();
            aligner = aligner.with_index(index.to_str().unwrap(), None).unwrap();
            Ok(aligner)
        }
    }

    fn __bool__(&self) -> PyResult<bool> {
        Ok(!self.idx.is_none())
    }

    #[getter]
    fn seq_names(&self) -> PyResult<Vec<String>> {
        if !self._has_index().unwrap() {
            return Err(PyRuntimeError::new_err("Index hasn't loaded"));
        }
        unsafe {
            let ns = (self.idx.unwrap()).n_seq;
            let mut sn = Vec::with_capacity(ns as usize);
            for i in 0..ns {
                sn.push(
                    std::ffi::CStr::from_ptr((*((self.idx.unwrap()).seq.offset(i as isize))).name)
                        .to_str()
                        .unwrap()
                        .to_string(),
                )
            }
            Ok(sn)
        }
    }

    /// Align a sequence with optional Metadata tuple. If provided, the tuple MUST be in the shape of (read_number: int, channel_number: int, read_id: str)
    fn _align(&self, seq: String, metadata: Option<&PyTuple>) -> PyResult<AlignmentResult> {
        let metadata = match metadata {
            Some(m) => {
                let inner: (i32, i32, String) = m.extract()?;
                MetaData::from(inner)
            }
            None => MetaData::from((0, 0, String::from(""))),
        };
        // do the heavy work
        let mappings = self.map(seq.as_bytes(), false, false, None, None).unwrap();

        // self._pool.join();
        // let alignment
        let return_metadata: (i32, i32, String) = (
            metadata.read_number,
            metadata.channel_number,
            metadata.read_id,
        );
        Ok(AlignmentResult {
            metadata: return_metadata,
            mappings: mappings.into_iter(),
        })
    }

    /// Align a sequence with optional Metadata tuple. If provided, the tuple MUST be in the shape of (read_number: int, channel_number: int, read_id: str)
    fn align_batch(
        &self,
        seqs: &PyIterator,
    ) -> PyResult<AlignmentBatchResultIter> {
        let mut res = AlignmentBatchResultIter::new(self.map_threads);
        let wq = Arc::clone(&self.work_queue);

        for seq_tup in seqs.iter()?
         {
            let (seq, r_n, c_n, r_id): (String, u32, u32, String) = seq_tup?.extract().unwrap();
            // println!("pushing rn {}", seq);
            wq.push(WorkQueue::Work(seq)).unwrap();
            res.sequences_sent += 1;
        }

        for _ in 0..self.map_threads {
            wq.push(WorkQueue::Done).unwrap();
        }
        // do the heavy work
        self.map_thread(&mut res)?;
        // let alignment
        // let return_metadata: (i32, i32, String) = (metadata.read_number, metadata.channel_number, String::from("hdea"));
        return Ok(res);
    }
}
/// This is for demonstrating how to return a value from __next__
#[pyclass]
pub struct AlignmentBatchResultIter {
    results_queue: Arc<ArrayQueue<WorkQueue<AlignmentResult>>>,
    finished_threads: u8,
    sequences_aligned: Arc<Mutex<u32>>,
    sequences_sent: u32,
    _n_threads: u8,
}

/// Iterator for the batch results from a multi threaded call to mapper
#[pymethods]
impl AlignmentBatchResultIter {
    #[new]
    pub fn new(n_threads: usize) -> Self {
        AlignmentBatchResultIter {
            results_queue: Arc::new(ArrayQueue::<WorkQueue<AlignmentResult>>::new(50000)),
            finished_threads: 0,
            sequences_aligned: Arc::new(Mutex::new(0)),
            sequences_sent: 0,
            _n_threads: n_threads as u8,
        }
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> IterNextOutput<AlignmentResult, &'static str> {
        let rq = Arc::clone(&self.results_queue);
        let mut work_item = rq.pop();
        // Not sure baout the next couple of lines of code
        let backoff = crossbeam::utils::Backoff::new();
        // let p = Parker::new();

        // println!("{:#?}", work_item);
        let mut mapped_work = None;
        loop {
            match work_item {
                Some(WorkQueue::Work(_)) => {
                    mapped_work = work_item;
                    break;
                }
                None => {
                    work_item = rq.pop();
                }
                Some(WorkQueue::Done) => {
                    self.finished_threads += 1;
                    work_item = rq.pop();
                }
                _ => {
                    work_item = rq.pop();
                }
            };
            // if *self.sequences_aligned.lock().unwrap() == self.sequences_sent && self.finished_threads == 8 {
            if self.finished_threads == self._n_threads {
                mapped_work = Some(WorkQueue::Done);
                break;
            }
        }

        return match mapped_work {
            Some(WorkQueue::Work(ar)) => IterNextOutput::Yield(ar),
            Some(WorkQueue::Done) => IterNextOutput::Return("Finished"),
            Some(WorkQueue::Starting) => IterNextOutput::Return("Error "),
            None => IterNextOutput::Return("Error"),
        };
    }
}

impl Aligner {
    pub fn map_thread(
        &self,
        res: &mut AlignmentBatchResultIter,
    ) -> Result<(), PyErr> {
        // TODO: Make threads Builder Argument
        // let metadata = match metadata {
        //     Some(m) => {
        //         let inner: (i32, i32, String) = m.extract()?;
        //         MetaData::from(inner)
        //     },
        //     None => {
        //         MetaData::from((0, 0, String::from("")))
        //     }
        // };

        // res.sequences_sent = seqs.len() as u32;
        res.sequences_sent = 100 as u32;

        // 8 threads
        for _ in 0..self.map_threads {
            let work_queue = Arc::clone(&self.work_queue);
            let results_queue = Arc::clone(&res.results_queue);
            let counter = Arc::clone(&res.sequences_aligned);
            // This could/should be done in the new thread, but wanting to test out the ability to move...
            // let mut aligner = aligner.clone();
            let aligner = self.clone();
            std::thread::spawn(move || loop {
                // let backoff = crossbeam::utils::Backoff::new();
                let work = work_queue.pop();
                match work {
                    Some(WorkQueue::Work(sequence)) => {

                        let mappings = aligner
                            .map(sequence.as_bytes(), false, false, None, None)
                            .expect("Unable to align");
                        // println!("MAppings {:#?}", mappings);
                        let ar = AlignmentResult {
                            metadata: (1, 2, String::from("3")),
                            mappings: mappings.into_iter(),
                        };
                        results_queue.push(WorkQueue::Work(ar)).unwrap();
                        // let mut num = counter.lock().unwrap();
                        // *num += 1;
                        // mem::drop(num);
                    }
                    Some(WorkQueue::Done) => {
                        results_queue.push(WorkQueue::Done).unwrap();
                        break;
                    }
                    None => std::thread::sleep(Duration::from_millis(100)),
                    _ => std::thread::sleep(Duration::from_millis(100)),
                }
            });
        }

        // This thread feeds all the sequences we have seen into the work queue to be processed above

        Ok(())
    }

    fn _has_index(&self) -> PyResult<bool> {
        Ok(self.idx.is_some())
    }

    /// Ergonomic function for Aligner. Just to see if people prefer this over the
    /// preset() function.
    /// ```
    /// # use mappy_rs::*;
    /// Aligner::builder().map_ont();
    /// ```
    pub fn map_ont(self) -> Self {
        self.preset(Preset::MapOnt)
    }

    /// Create an aligner using a preset.
    pub fn preset(self, preset: Preset) -> Self {
        let mut idxopt = IdxOpt::default();
        let mut mapopt = MapOpt::default();

        unsafe {
            // Set preset
            mm_set_opt(preset.into(), &mut idxopt, &mut mapopt)
        };

        Self {
            idxopt,
            mapopt,
            ..Default::default()
        }
    }

    /// XOR each extra flag provided at initialisation of the aligner onto the map_opt 
    /// This should be called after present has been called - TODO implement check for this
    pub fn with_extra_flags(mut self, extra_flags: Vec<u64>) -> Self {
        for flag in extra_flags {
            self.mapopt.flag |= flag as i64;
        }
        self
    }

    /// Set Alignment mode / cigar mode in minimap2
    pub fn with_cigar(mut self) -> Self {
        // Make sure MM_F_CIGAR flag isn't already set
        assert!((self.mapopt.flag & MM_F_CIGAR as i64) == 0);

        self.mapopt.flag |= MM_F_CIGAR as i64;
        self
    }

    /// Set the number of mapping threads (prefer to use the struct config)
    pub fn with_map_threads(mut self, threads: usize) -> Self {
        self.map_threads = threads;
        self
    }

    /// Set the number of threads (prefer to use the struct config)
    pub fn with_idx_threads(mut self, threads: usize) -> Self {
        self.idx_threads = threads;
        self
    }
    /// Set index parameters for minimap2 using builder pattern
    /// Creates the index as well with the given number of threads (set at struct creation)
    ///
    /// Parameters:
    /// path: Location of pre-built index or FASTA/FASTQ file (may be gzipped or plaintext)
    /// Output: Option (None) or a filename
    ///
    /// Returns the aligner with the index set
    ///
    /// ```
    /// # use mappy_rs::*;
    /// // Do not save the index file
    /// Aligner::builder().map_ont().with_index("test_data/test_data.fasta", None);
    ///
    /// // Save the index file as my_index.mmi
    /// Aligner::builder().map_ont().with_index("test_data/test_data.fasta", Some("my_index.mmi"));
    ///
    /// // Use the previously built index
    /// Aligner::builder().map_ont().with_index("my_index.mmi", None);
    /// ```
    pub fn with_index(mut self, path: &str, output: Option<&str>) -> Result<Self, &'static str> {
        // Confirm file exists
        if !Path::new(path).exists() {
            return Err("File does not exist");
        }

        // Confirm file is not empty
        if Path::new(path).metadata().unwrap().len() == 0 {
            return Err("File is empty");
        }

        let path = match std::ffi::CString::new(path) {
            Ok(path) => path,
            Err(_) => return Err("Invalid path"),
        };

        let output = match output {
            Some(output) => match std::ffi::CString::new(output) {
                Ok(output) => output,
                Err(_) => return Err("Invalid output"),
            },
            None => std::ffi::CString::new(Vec::new()).unwrap(),
        };

        let idx_reader = MaybeUninit::new(unsafe {
            mm_idx_reader_open(path.as_ptr(), &self.idxopt, output.as_ptr())
        });

        let mut idx: MaybeUninit<*mut mm_idx_t> = MaybeUninit::uninit();

        let mut idx_reader = unsafe { idx_reader.assume_init() };

        unsafe {
            // Just a test read? Just following: https://github.com/lh3/minimap2/blob/master/python/mappy.pyx#L147
            idx = MaybeUninit::new(mm_idx_reader_read(
                // self.idx_reader.as_mut().unwrap() as *mut mm_idx_reader_t,
                &mut *idx_reader as *mut mm_idx_reader_t,
                self.idx_threads as libc::c_int,
            ));
            // Close the reader
            mm_idx_reader_close(idx_reader);
            // Set index opts
            mm_mapopt_update(&mut self.mapopt, *idx.as_ptr());
            // Idx index name
            mm_idx_index_name(idx.assume_init());
        }

        self.idx = Some(unsafe { *idx.assume_init() });
        Ok(self)
    }

    // https://github.com/lh3/minimap2/blob/master/python/mappy.pyx#L164
    // TODO: I doubt extra_flags is working properly...
    // TODO: Python allows for paired-end mapping with seq2: Option<&[u8]>, but more work to implement
    /// Aligns a given sequence (as bytes) to the index associated with this aligner
    ///
    /// Parameters:
    /// seq: Sequence to align
    /// cs: Whether to output CIGAR string
    /// MD: Whether to output MD tag
    /// max_frag_len: Maximum fragment length
    ///
    // Make sure index is set
    pub fn map(
        &self,
        seq: &[u8],
        cs: bool,
        md: bool, // TODO
        max_frag_len: Option<usize>,
        extra_flags: Option<Vec<u64>>,
    ) -> Result<Vec<Mapping>, &'static str> {
        // Make sure index is set
        if !self.has_index() {
            return Err("No index");
        }

        // Make sure sequence is not empty
        if seq.len() == 0 {
            return Err("Sequence is empty");
        }

        let mut mm_reg: MaybeUninit<*mut mm_reg1_t> = MaybeUninit::uninit();

        // Number of results
        let mut n_regs: i32 = 0;
        let mut map_opt = self.mapopt.clone();

        // if max_frag_len is not None: map_opt.max_frag_len = max_frag_len
        if let Some(max_frag_len) = max_frag_len {
            map_opt.max_frag_len = max_frag_len as i32;
        }

        // if extra_flags is not None: map_opt.flag |= extra_flags
        if let Some(extra_flags) = extra_flags {
            for flag in extra_flags {
                map_opt.flag |= flag as i64;
            }
        }

        let mappings = BUF.with(|buf| {
            let km = unsafe { mm_tbuf_get_km(buf.borrow_mut().buf) };

            mm_reg = MaybeUninit::new(unsafe {
                mm_map(
                    self.idx.as_ref().unwrap() as *const mm_idx_t,
                    seq.len() as i32,
                    seq.as_ptr() as *const i8,
                    &mut n_regs,
                    buf.borrow_mut().buf,
                    &mut map_opt,
                    std::ptr::null(),
                )
            });

            let mut mappings = Vec::with_capacity(n_regs as usize);

            for i in 0..n_regs {
                unsafe {
                    let reg_ptr = (*mm_reg.as_ptr()).offset(i as isize);
                    // println!("{:#?}", *reg_ptr);
                    let const_ptr = reg_ptr as *const mm_reg1_t;
                    let reg: mm_reg1_t = *reg_ptr;

                    // TODO: Get all contig names and store as Cow<String> somewhere centralized...
                    let contig: *mut ::std::os::raw::c_char =
                        (*(self.idx.unwrap()).seq.offset(reg.rid as isize)).name;

                    let alignment = if !reg.p.is_null() {
                        let p = &*reg.p;
                        let n_cigar = p.n_cigar;
                        let cigar: Vec<u32> = p.cigar.as_slice(n_cigar as usize).to_vec();
                        if cs {
                            // let mut cs_string: *mut std::ffi::c_char = std::ptr::null_mut();
                            let mut cs_string: *mut libc::c_char = std::ptr::null_mut();
                            let mut m_cs_string: libc::c_int = 0i32;

                            let cs_len = mm_gen_cs(
                                km,
                                &mut cs_string,
                                &mut m_cs_string,
                                &self.idx.unwrap() as *const mm_idx_t,
                                const_ptr,
                                seq.as_ptr() as *const i8,
                                true.into(),
                            );

                            let cs_string = std::ffi::CStr::from_ptr(cs_string)
                                .to_str()
                                .unwrap()
                                .to_string();

                            Some(Alignment {
                                is_primary: true,
                                cigar: Some(format!("cs:Z::{}", cs_string)),
                            })
                        } else {
                            Some(Alignment {
                                is_primary: true,
                                cigar: None,
                            })
                        }
                    } else {
                        None
                    };

                    mappings.push(Mapping {
                        target_name: Some(
                            std::ffi::CStr::from_ptr(contig)
                                .to_str()
                                .unwrap()
                                .to_string(),
                        ),
                        target_len: (*(self.idx.unwrap()).seq.offset(reg.rid as isize)).len as i32,
                        target_start: reg.rs,
                        target_end: reg.re,
                        query_name: None,
                        query_len: Some(seq.len() as u32),
                        query_start: reg.qs,
                        query_end: reg.qe,
                        strand: if reg.rev() == 0 {
                            Strand::Forward
                        } else {
                            Strand::Reverse
                        },
                        match_len: reg.mlen,
                        block_len: reg.blen,
                        mapq: reg.mapq(),
                        alignment,
                    });
                }
            }

            mappings
        });
        Ok(mappings)
    }

    // This is in the python module, so copied here...
    pub fn has_index(&self) -> bool {
        self.idx.is_some()
    }
}

// impl Drop for Aligner {
//     fn drop(&mut self) {
//         if self.idx.is_some() {
//             unsafe { mm_idx_destroy(self.idx.unwrap()) };
//         }
//     }
// }

#[derive(PartialEq, Eq)]
pub enum FileFormat {
    FASTA,
    FASTQ,
}

#[allow(dead_code)]
pub fn detect_file_format(buffer: &[u8]) -> Result<FileFormat, &'static str> {
    let buffer = from_utf8(&buffer).expect("Unable to parse file as UTF-8");
    if buffer.starts_with(">") {
        Ok(FileFormat::FASTA)
    } else if buffer.starts_with("@") {
        Ok(FileFormat::FASTQ)
    } else {
        Err("Unknown file format")
    }
}

#[derive(PartialEq, Eq, Debug)]
pub enum CompressionType {
    GZIP,
    BZIP2,
    XZ,
    RAR,
    ZSTD,
    LZ4,
    LZMA,
    NONE,
}

/// Return the compression type of a file
#[allow(dead_code)]
pub fn detect_compression_format(buffer: &[u8]) -> Result<CompressionType, &'static str> {
    Ok(match buffer {
        [0x1F, 0x8B, ..] => CompressionType::GZIP,
        [0x42, 0x5A, ..] => CompressionType::BZIP2,
        [0xFD, b'7', b'z', b'X', b'Z', 0x00] => CompressionType::XZ,
        [0x28, 0xB5, 0x2F, 0xFD, ..] => CompressionType::LZMA,
        [0x5D, 0x00, ..] => CompressionType::LZMA,
        [0x1F, 0x9D, ..] => CompressionType::LZMA,
        [0x37, 0x7A, 0xBC, 0xAF, 0x27, 0x1C] => CompressionType::ZSTD,
        [0x04, 0x22, 0x4D, 0x18, ..] => CompressionType::LZ4,
        [0x08, 0x22, 0x4D, 0x18, ..] => CompressionType::LZ4,
        [0x52, 0x61, 0x72, 0x21, 0x1A, 0x07] => CompressionType::RAR,
        _ => CompressionType::NONE,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::mem::MaybeUninit;

    #[test]
    fn does_it_work() {
        let mut mm_idxopt = MaybeUninit::uninit();
        let mut mm_mapopt = MaybeUninit::uninit();

        unsafe { mm_set_opt(&0, mm_idxopt.as_mut_ptr(), mm_mapopt.as_mut_ptr()) };
    }

    #[test]
    fn create_index_file_missing() {
        let result = Aligner::builder()
            .preset(Preset::MapOnt)
            .with_idx_threads(2)
            .with_index(
                "test_data/test.fa_FILE_NOT_FOUND",
                Some("test_FILE_NOT_FOUND.mmi"),
            );
        assert!(result.is_err());
    }

    // #[test]
    // fn create_index() {
    //     let mut aligner = Aligner::preset(Preset::MapOnt).with_threads(1);

    //     println!("{}", aligner.idxopt.w);

    //     assert!(aligner.idxopt.w == 10);

    //     aligner = aligner
    //         .with_index("test_data/test_data.fasta", Some("test.mmi"))
    //         .unwrap();
    // }

    #[test]
    fn test_mapping() {
        let mut aligner = Aligner::builder()
            .preset(Preset::MapOnt)
            .with_idx_threads(2);

        aligner = aligner.with_index("tests/yeast_ref.mmi", None).unwrap();

        aligner
            .map(
                "ACGGTAGAGAGGAAGAAGAAGGAATAGCGGACTTGTGTATTTTATCGTCATTCGTGGTTATCATATAGTTTATTGATTTGAAGACTACGTAAGTAATTTGAGGACTGATTAAAATTTTCTTTTTTAGCTTAGAGTCAATTAAAGAGGGCAAAATTTTCTCAAAAGACCATGGTGCATATGACGATAGCTTTAGTAGTATGGATTGGGCTCTTCTTTCATGGATGTTATTCAGAAGGAGTGATATATCGAGGTGTTTGAAACACCAGCGACACCAGAAGGCTGTGGATGTTAAATCGTAGAACCTATAGACGAGTTCTAAAATATACTTTGGGGTTTTCAGCGATGCAAAA".as_bytes(),
                false,
                false,
                None,
                None,
            )
            .unwrap();
        let mappings = aligner.map("ACGGTAGAGAGGAAGAAGAAGGAATAGCGGACTTGTGTATTTTATCGTCATTCGTGGTTATCATATAGTTTATTGATTTGAAGACTACGTAAGTAATTTGAGGACTGATTAAAATTTTCTTTTTTAGCTTAGAGTCAATTAAAGAGGGCAAAATTTTCTCAAAAGACCATGGTGCATATGACGATAGCTTTAGTAGTATGGATTGGGCTCTTCTTTCATGGATGTTATTCAGAAGGAGTGATATATCGAGGTGTTTGAAACACCAGCGACACCAGAAGGCTGTGGATGTTAAATCGTAGAACCTATAGACGAGTTCTAAAATATACTTTGGGGTTTTCAGCGATGCAAAA".as_bytes(), false, false, None, None).unwrap();
        println!("{:#?}", mappings);

        let aligner = aligner.with_cigar();

        aligner
            .map(
                "ATGAGCAAAATATTCTAAAGTGGAAACGGCACTAAGGTGAACTAAGCAACTTAGTGCAAAAc".as_bytes(),
                true,
                false,
                None,
                None,
            )
            .unwrap();

        let mappings = aligner.map("atCCTACACTGCATAAACTATTTTGcaccataaaaaaaagttatgtgtgGGTCTAAAATAATTTGCTGAGCAATTAATGATTTCTAAATGATGCTAAAGTGAACCATTGTAatgttatatgaaaaataaatacacaattaagATCAACACAGTGAAATAACATTGATTGGGTGATTTCAAATGGGGTCTATctgaataatgttttatttaacagtaatttttatttctatcaatttttagtaatatctacaaatattttgttttaggcTGCCAGAAGATCGGCGGTGCAAGGTCAGAGGTGAGATGTTAGGTGGTTCCACCAACTGCACGGAAGAGCTGCCCTCTGTCATTCAAAATTTGACAGGTACAAACAGactatattaaataagaaaaacaaactttttaaaggCTTGACCATTAGTGAATAGGTTATATGCTTATTATTTCCATTTAGCTTTTTGAGACTAGTATGATTAGACAAATCTGCTTAGttcattttcatataatattgaGGAACAAAATTTGTGAGATTTTGCTAAAATAACTTGCTTTGCTTGTTTATAGAGGCacagtaaatcttttttattattattataattttagattttttaatttttaaat".as_bytes(), true, false, None, None).unwrap();
        println!("{:#?}", mappings);
    }

    #[test]
    fn test_aligner_config_and_mapping() {
        let mut aligner = Aligner::builder().map_ont().with_idx_threads(2);
        aligner = aligner
            .with_index("tests/yeast_ref.mmi", Some("test.mmi"))
            .unwrap()
            .with_cigar();

        aligner
            .map(
                "ATGAGCAAAATATTCTAAAGTGGAAACGGCACTAAGGTGAACTAAGCAACTTAGTGCAAAAc".as_bytes(),
                true,
                false,
                None,
                None,
            )
            .unwrap();
        let mappings = aligner.map("atCCTACACTGCATAAACTATTTTGcaccataaaaaaaagGGACatgtgtgGGTCTAAAATAATTTGCTGAGCAATTAATGATTTCTAAATGATGCTAAAGTGAACCATTGTAatgttatatgaaaaataaatacacaattaagATCAACACAGTGAAATAACATTGATTGGGTGATTTCAAATGGGGTCTATctgaataatgttttatttaacagtaatttttatttctatcaatttttagtaatatctacaaatattttgttttaggcTGCCAGAAGATCGGCGGTGCAAGGTCAGAGGTGAGATGTTAGGTGGTTCCACCAACTGCACGGAAGAGCTGCCCTCTGTCATTCAAAATTTGACAGGTACAAACAGactatattaaataagaaaaacaaactttttaaaggCTTGACCATTAGTGAATAGGTTATATGCTTATTATTTCCATTTAGCTTTTTGAGACTAGTATGATTAGACAAATCTGCTTAGttcattttcatataatattgaGGAACAAAATTTGTGAGATTTTGCTAAAATAACTTGCTTTGCTTGTTTATAGAGGCacagtaaatcttttttattattattataattttagattttttaatttttaaat".as_bytes(), true, true, None, None).unwrap();
        println!("{:#?}", mappings);
    }
}

#[pymodule]
fn mappy_rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Aligner>()?;
    m.add_class::<Mapping>()?;
    m.add_class::<AlignmentResult>()?;
    m.add_class::<AlignmentBatchResultIter>()?;
    Ok(())
}
