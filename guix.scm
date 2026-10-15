;; Guix package for the fast_trimming pipeline.
;;   Build:            guix build -f guix.scm
;;   Isolated shell:   guix shell --container -f guix.scm
;;   Pinned:           guix time-machine -C channels.scm -- shell --container -f guix.scm
(use-modules (guix packages)
             (guix gexp)
             (guix build-system gnu)
             ((guix licenses) #:prefix license:)
             (gnu packages base)
             (gnu packages bash)
             (gnu packages bioinformatics)
             (gnu packages compression)
             (gnu packages gawk)
             (gnu packages linux)
             (gnu packages web))

(define %source-dir (dirname (current-filename)))

;; Only take the scripts, src/ and the test data, never other
;; sequencing data or outputs lying next to them.
(define (source-file? file stat)
  (let ((name (basename file)))
    (or (eq? (stat:type stat) 'directory)
        (string-suffix? ".c" name)
        (string-suffix? ".h" name)
        (string-suffix? ".sh" name)
        (string=? "makefile" name)
        (string=? "test-data.fastq.gz" name)
        (string=? "test-data.counts.orig" name))))

(package
  (name "fast-trimming")
  (version "0.1")
  (source (local-file %source-dir "fast-trimming-checkout"
                      #:recursive? #t
                      #:select? source-file?))
  (build-system gnu-build-system)
  (arguments
   (list
    #:make-flags #~(list "CC=gcc" (string-append "PREFIX=" #$output))
    #:phases
    #~(modify-phases %standard-phases
        (delete 'configure)
        (add-after 'unpack 'chdir
          (lambda _ (chdir "src")))
        (add-after 'install 'wrap-scripts
          ;; Makes the scripts find their tools in the store only,
          ;; whatever the caller's PATH contains.
          (lambda* (#:key inputs #:allow-other-keys)
            (let ((bin (string-append #$output "/bin")))
              (for-each
               (lambda (script)
                 (wrap-program (string-append bin "/" script)
                   `("PATH" ":" prefix
                     ,(cons bin
                            (map (lambda (cmd)
                                   (dirname (search-input-file
                                             inputs (string-append "bin/" cmd))))
                                 '("cutadapt" "vsearch" "jq" "awk" "rev"
                                   "zcat" "tee" "wc" "date"))))))
               '("fast_trimming.sh" "comp_dna.sh")))))
        ;; The pipeline can only be tested once installed, wrapped and
        ;; its shebangs patched, so the test runs after install instead
        ;; of before it.
        (delete 'check)
        (add-after 'patch-shebangs 'check
          ;; Runs the installed pipeline on the test data and fails the
          ;; build if its counts differ from the reference ones.
          (lambda* (#:key tests? #:allow-other-keys)
            (when tests?
              (let ((reference (canonicalize-path "../test-data.counts.orig")))
                (mkdir "test-run")
                (copy-file "../test-data.fastq.gz" "test-run/test-data.fastq.gz")
                (with-directory-excursion "test-run"
                  (invoke (string-append #$output "/bin/fast_trimming.sh")
                          "-t" "0.1" "-is" "33"
                          "-a1" "CCACCAGAGT" "-a2" "GCGCAGACCG"
                          "-i" "test-data.fastq.gz")
                  (invoke "diff" "test-data.counts" reference)))))))))
  (inputs (list bash-minimal coreutils cutadapt gawk gzip jq util-linux
                vsearch zlib))
  (home-page "https://github.com/timothevanmeter/trimming")
  (synopsis "Trim NGS reads to a target insert and count unique variants")
  (description "Extracts the sequence between two adapters (constant sequences flanking the region of interest), filters it by length, quality and stop codons, and counts each unique amino acid variant.")
  (license license:gpl3+))
