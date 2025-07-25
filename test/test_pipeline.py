#!/usr/bin/env python3


"""
Test module for pipeline.py
"""

import unittest
from unittest.mock import ANY as mock_ANY
from unittest.mock import patch, call, MagicMock, mock_open, Mock
from scripts.pipeline import *



class TestBwaAlign(unittest.TestCase):
    """
    Unit tests for the bwa_align function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_mkdir = patch('os.mkdir').start()
        self.mock_system = patch('os.system').start()
        self.mock_log_code = patch('scripts.pipeline.log_code').start()

        self.sra_id = "test_sra_id_bwa"
        self.ref_genome = "/data/ref/genome.fasta"
        self.output_dir = "/tmp/bwa_output"
        self.log_scr = "/var/log/bwa_script.log"
        self.threads = 4
        self.expected_sam = os.path.join(
            self.output_dir,
            self.sra_id + '.sam'
        )

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_bwa_align_single_fastq(self):
        l_fastq = ["/data/reads/read1.fastq.gz"]
        result = bwa_align(
            self.sra_id,
            l_fastq,
            self.ref_genome,
            self.output_dir,
            self.log_scr,
            self.threads
        )

        self.mock_mkdir.assert_called_once_with(self.output_dir)
        expected_cmd = (
            f'bwa mem -M -t {self.threads} {self.ref_genome}' +\
                f' {l_fastq[0]} > {self.expected_sam}'
        )
        self.mock_log_code.assert_called_once_with(
            expected_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_called_once_with(expected_cmd)
        self.assertEqual(result, self.expected_sam)

    def test_successful_bwa_align_paired_fastq(self):
        l_fastq = [
            "/data/reads/read1.fastq.gz",
            "/data/reads/read2.fastq.gz"
        ]
        result = bwa_align(
            self.sra_id,
            l_fastq,
            self.ref_genome,
            self.output_dir,
            self.log_scr,
            self.threads
        )

        self.mock_mkdir.assert_called_once_with(self.output_dir)
        expected_cmd = (
            f'bwa mem -M -t {self.threads} {self.ref_genome}' +\
                f' {l_fastq[0]} {l_fastq[1]} > {self.expected_sam}'
        )
        self.mock_log_code.assert_called_once_with(
            expected_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_called_once_with(expected_cmd)
        self.assertEqual(result, self.expected_sam)

    def test_bwa_align_directory_already_exists(self):
        self.mock_mkdir.side_effect = FileExistsError
        l_fastq = ["/data/reads/read1.fastq.gz"]
        result = bwa_align(
            self.sra_id,
            l_fastq,
            self.ref_genome,
            self.output_dir,
            self.log_scr,
            self.threads
        )

        self.mock_mkdir.assert_called_once_with(self.output_dir)
        expected_cmd = (
            f'bwa mem -M -t {self.threads} {self.ref_genome}' +\
                f' {l_fastq[0]} > {self.expected_sam}'
        )
        self.mock_log_code.assert_called_once_with(
            expected_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_called_once_with(expected_cmd)
        self.assertEqual(result, self.expected_sam)


class TestChangeOutputOwnership(unittest.TestCase):
    """
    Unit tests for the change_output_ownership function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_os_environ = patch(
            'os.environ',
            new_callable=MagicMock
        ).start()
        self.mock_os_walk = patch('os.walk').start()
        self.mock_os_chown = patch('os.chown').start()
        self.mock_log_print = patch('scripts.pipeline.log_print').start()

        self.output_dir = "/test/output/directory"
        self.log_file = "/var/log/ownership.log"

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_ownership_change(self):
        """
        Tests the successful scenario where ownership is changed for all 
        files and directories.
        """
        # Set up environment variables for the test
        self.mock_os_environ.__getitem__.side_effect = lambda key: {
            "HOST_UID": "1000", "HOST_GID": "1000"
        }[key]

        # Simulate directory structure
        self.mock_os_walk.return_value = [
            (
                self.output_dir,
                ['subdir1', 'subdir2'],
                ['file1.txt', 'file2.log']
            ),
            (
                os.path.join(self.output_dir, 'subdir1'),
                [],
                ['subfile1.txt']
            ),
            (
                os.path.join(self.output_dir, 'subdir2'),
                [],
                []
            )
        ]

        change_output_ownership(self.output_dir, self.log_file)

        # Assert os.chown calls
        expected_chown_calls = [
            call(os.path.join(self.output_dir, 'subdir1'), 1000, 1000),
            call(os.path.join(self.output_dir, 'subdir2'), 1000, 1000),
            call(os.path.join(self.output_dir, 'file1.txt'), 1000, 1000),
            call(os.path.join(self.output_dir, 'file2.log'), 1000, 1000),
            call(
                os.path.join(self.output_dir, 'subdir1', 'subfile1.txt'),
                1000,
                1000
            )
        ]
        self.mock_os_chown.assert_has_calls(
            expected_chown_calls,
            any_order=True
        )
        self.assertEqual(
            self.mock_os_chown.call_count,
            len(expected_chown_calls)
        )

        # Assert log_print calls for successful changes
        expected_log_prints = [
            call(
                f"Changed ownership of file " +\
                    f"{os.path.join(self.output_dir, 'subdir1')}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f"Changed ownership of file " +\
                    f"{os.path.join(self.output_dir, 'subdir2')}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f"Changed ownership of file " +\
                    f"{os.path.join(self.output_dir, 'file1.txt')}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f"Changed ownership of file " +\
                    f"{os.path.join(self.output_dir, 'file2.log')}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f"Changed ownership of file " +\
                    os.path.join(
                        self.output_dir,
                        'subdir1',
                        'subfile1.txt'
                    ),
                level='info',
                log_file=self.log_file
            )
        ]
        self.mock_log_print.assert_has_calls(
            expected_log_prints,
            any_order=True
        )
        self.assertEqual(
            self.mock_log_print.call_count,
            len(expected_log_prints)
        )

    def test_key_error_no_env_vars(self):
        """
        Tests the scenario where HOST_UID or HOST_GID environment 
        variables are not set.
        Ensures a warning is logged and chown is not called.
        """
        # Configure os.environ to raise KeyError when accessing HOST_UID
        self.mock_os_environ.__getitem__.side_effect = KeyError("HOST_UID")

        change_output_ownership(self.output_dir, self.log_file)

        # Assert that os.chown was never called
        self.mock_os_chown.assert_not_called()

        # Assert that a warning message is logged
        expected_warning_message = (
            'Permissions cannot be changed.\nTo do it manually, run:'
            f'\nsudo chmod 777 {self.output_dir}'
        )
        self.mock_log_print.assert_called_once_with(
            expected_warning_message,
            level='warn',
            log_file=self.log_file
        )
        self.assertEqual(self.mock_log_print.call_count, 1)

    def test_permission_error_during_chown(self):
        """
        Tests the scenario where os.chown raises a PermissionError for 
        some files.
        Ensures a warning is logged for the specific file and execution 
        continues.
        """
        self.mock_os_environ.__getitem__.side_effect = lambda key: {
            "HOST_UID": "1000",
            "HOST_GID": "1000"
        }[key]

        # Simulate directory structure with one file causing 
        # PermissionError
        file_with_permission_error = os.path.join(
            self.output_dir,
            'restricted_file.txt'
        )
        self.mock_os_walk.return_value = [
            (self.output_dir, [], ['file1.txt', 'restricted_file.txt'])
        ]

        # Configure os.chown to raise PermissionError for a specific file
        def chown_side_effect(path, uid, gid):
            if path == file_with_permission_error:
                raise PermissionError(
                    f"Operation not permitted: '{path}'"
                )
            # For other files, do nothing (simulate success)
            pass

        self.mock_os_chown.side_effect = chown_side_effect

        change_output_ownership(self.output_dir, self.log_file)

        # Assert chown was attempted for both files
        self.assertEqual(self.mock_os_chown.call_count, 2)
        self.mock_os_chown.assert_any_call(
            os.path.join(self.output_dir, 'file1.txt'),
            1000,
            1000
        )
        self.mock_os_chown.assert_any_call(
            file_with_permission_error,
            1000,
            1000
        )

        # Assert log_print calls: one for success, one for warning
        self.assertEqual(self.mock_log_print.call_count, 2)
        self.mock_log_print.assert_any_call(
            "Changed ownership of file " +\
                f"{os.path.join(self.output_dir, 'file1.txt')}",
            level='info',
            log_file=self.log_file
        )
        self.mock_log_print.assert_any_call(
            f"Permission denied on {file_with_permission_error}: " +\
                f"Operation not permitted: '{file_with_permission_error}'",
            level='warn',
            log_file=self.log_file
        )

    def test_file_not_found_error_during_chown(self):
        """
        Tests the scenario where os.chown raises a FileNotFoundError for 
        some files.
        Ensures the error is caught and execution continues without 
        logging a warning.
        """
        self.mock_os_environ.__getitem__.side_effect = lambda key: {
            "HOST_UID": "1000",
            "HOST_GID": "1000"
        }[key]

        # Simulate directory structure with one file disappearing
        file_not_found = os.path.join(self.output_dir, 'missing_file.txt')
        self.mock_os_walk.return_value = [
            (self.output_dir, [], ['file1.txt', 'missing_file.txt'])
        ]

        # Configure os.chown to raise FileNotFoundError for a specific file
        def chown_side_effect(path, uid, gid):
            if path == file_not_found:
                raise FileNotFoundError(
                    f"No such file or directory: '{path}'"
                )
            pass

        self.mock_os_chown.side_effect = chown_side_effect

        change_output_ownership(self.output_dir, self.log_file)

        # Assert chown was attempted for both files
        self.assertEqual(self.mock_os_chown.call_count, 2)
        self.mock_os_chown.assert_any_call(
            os.path.join(self.output_dir, 'file1.txt'),
            1000,
            1000
        )
        self.mock_os_chown.assert_any_call(file_not_found, 1000, 1000)

        # Assert log_print calls: only for the successful change, not for 
        # FileNotFoundError
        self.assertEqual(self.mock_log_print.call_count, 1)
        self.mock_log_print.assert_called_once_with(
            f"Changed ownership of file " +\
                f"{os.path.join(self.output_dir, 'file1.txt')}",
            level='info',
            log_file=self.log_file
        )


class TestCompressIndexVcf(unittest.TestCase):
    """
    Unit tests for the compress_index_vcf function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_mkdir = patch('os.mkdir').start()
        self.mock_system = patch('os.system').start()
        self.mock_log_print = patch('scripts.pipeline.log_print').start()
        self.mock_log_code = patch('scripts.pipeline.log_code').start()

        self.sra_id = "test_sra_id_vcf"
        self.vcf_file = "/data/variants.vcf"
        self.output_dir = "/tmp/vcf_output"
        self.log_file = "/var/log/vcf_pipeline.log"
        self.log_scr = "/var/log/vcf_script.log"
        self.expected_compressed_vcf = os.path.join(
            self.output_dir,
            self.sra_id + '.vcf.gz'
        )
        self.expected_tbi_file = self.expected_compressed_vcf + '.tbi'

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_compression_and_indexing(self):
        # Simulate success for both bgzip and tabix
        self.mock_system.return_value = 0

        result = compress_index_vcf(
            self.sra_id,
            self.vcf_file,
            self.output_dir,
            self.log_file,
            self.log_scr
        )

        self.mock_mkdir.assert_called_once_with(self.output_dir)

        # Assert bgzip commands and logs
        expected_bgzip_cmd = f'bgzip -c {self.vcf_file} > ' +\
            f'{self.expected_compressed_vcf}'
        self.mock_log_print.assert_any_call(
            f"Compressing {self.vcf_file} with bgzip...",
            level='info',
            log_file=self.log_file
        )
        self.mock_log_code.assert_any_call(
            expected_bgzip_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_any_call(expected_bgzip_cmd)
        self.mock_log_print.assert_any_call(
            f"Compressed VCF: {self.expected_compressed_vcf}",
            level='info',
            log_file=self.log_file
        )

        # Assert tabix commands and logs
        expected_tabix_cmd = f'tabix -p vcf {self.expected_compressed_vcf}'
        self.mock_log_print.assert_any_call(
            f"Indexing {self.expected_compressed_vcf} with tabix...",
            level='info',
            log_file=self.log_file
        )
        self.mock_log_code.assert_any_call(
            expected_tabix_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_any_call(expected_tabix_cmd)
        self.mock_log_print.assert_any_call(
            f"VCF index: {self.expected_tbi_file}",
            level='info',
            log_file=self.log_file
        )

        self.assertEqual(result, self.expected_compressed_vcf)
        self.assertEqual(self.mock_system.call_count, 2)
        self.assertEqual(self.mock_log_print.call_count, 4)
        self.assertEqual(self.mock_log_code.call_count, 2)

    def test_compress_index_vcf_directory_already_exists(self):
        self.mock_mkdir.side_effect = FileExistsError
        self.mock_system.return_value = 0

        compress_index_vcf(
            self.sra_id,
            self.vcf_file,
            self.output_dir,
            self.log_file,
            self.log_scr
        )
        self.mock_mkdir.assert_called_once_with(self.output_dir)
        # Verify that the rest of the process still runs
        self.assertEqual(self.mock_system.call_count, 2)

    def test_tabix_command_failure(self):
        # bgzip succeeds, tabix fails
        self.mock_system.side_effect = [0, 1]

        result = compress_index_vcf(
            self.sra_id,
            self.vcf_file,
            self.output_dir,
            self.log_file,
            self.log_scr
        )

        expected_bgzip_cmd = f'bgzip -c {self.vcf_file} > ' +\
            f'{self.expected_compressed_vcf}'
        expected_tabix_cmd = f'tabix -p vcf {self.expected_compressed_vcf}'
        self.mock_system.assert_has_calls([
            call(expected_bgzip_cmd),
            call(expected_tabix_cmd)
        ])
        self.assertEqual(self.mock_system.call_count, 2)
        # All logs should still happen
        self.assertEqual(self.mock_log_print.call_count, 4)
        self.assertEqual(self.mock_log_code.call_count, 2)
        # Still returns path even if index not created
        self.assertEqual(result, self.expected_compressed_vcf)


class TestCreateCounts(unittest.TestCase):
    """
    Unit tests for the create_counts function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_open = patch('builtins.open', mock_open()).start()
        self.mock_os_environ = patch(
            'os.environ',
            new_callable=MagicMock
        ).start()
        self.mock_os_system = patch('os.system').start()

        # Define common test parameters
        self.counts_file = "/tmp/test_counts.txt"
        self.vcf_compressed = "/data/variants.vcf.gz"

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_counts_creation_multiple_regions(self):
        """
        Tests the creation of a counts file with multiple regions.
        Verifies file writes, environment variable settings, and system 
        calls.
        """
        regions = [
            ["chr1", 100, 200, "geneA"],
            ["chr2", 300, 400, "geneB"],
            ["chrX", 500, 600, "geneC"]
        ]

        # Simulate os.system's output being appended to the file

        create_counts(self.counts_file, regions, self.vcf_compressed)

        # Verify 'open' calls: first 'w' then 'a+' for each region
        self.mock_open.assert_has_calls([
            call(self.counts_file, "w"),
            call().__enter__(),
            call().write(""),
            call().__exit__(None, None, None),
            # For each region:
            call(self.counts_file, "a+"),
            call().__enter__(),
            call().write("geneA\t"),
            call().__exit__(None, None, None),
            call(self.counts_file, "a+"),
            call().__enter__(),
            call().write("geneB\t"),
            call().__exit__(None, None, None),
            call(self.counts_file, "a+"),
            call().__enter__(),
            call().write("geneC\t"),
            call().__exit__(None, None, None),
        ], any_order=False)

        # Verify os.environ settings for each region
        expected_environ_calls = []
        for chrom, start, end, name in regions:
            region_str = f"{chrom}:{start}-{end}"
            expected_environ_calls.extend([
                call.__setitem__('region_str', region_str),
                call.__setitem__('vcf_file', self.vcf_compressed),
                call.__setitem__('tmp_output', self.counts_file),
            ])
        self.mock_os_environ.assert_has_calls(
            expected_environ_calls,
            any_order=False
        )
        self.assertEqual(
            self.mock_os_environ.__setitem__.call_count,
            len(regions) * 3
        )

        # Verify os.system calls for each region
        expected_system_calls = []
        for chrom, start, end, name in regions:
            region_str = f"{chrom}:{start}-{end}"
            cmd = (
                'bcftools view -r "$region_str"'
                ' "$vcf_file" | grep -vc "^#" >> "$tmp_output"'
            )
            expected_system_calls.append(call(cmd))

        self.mock_os_system.assert_has_calls(
            expected_system_calls,
            any_order=False
        )
        self.assertEqual(self.mock_os_system.call_count, len(regions))

    def test_empty_regions_list(self):
        """
        Tests the scenario where the regions list is empty.
        Ensures the file is cleared and no further operations occur.
        """
        regions = []

        create_counts(self.counts_file, regions, self.vcf_compressed)

        # Verify 'open' is called only once to clear the file
        self.mock_open.assert_has_calls([
            call(self.counts_file, "w"),
            call().__enter__(),
            call().write(""),
            call().__exit__(None, None, None),
        ])
        self.assertEqual(self.mock_open.call_count, 1) 

        # Verify no os.environ or os.system calls
        self.mock_os_environ.assert_not_called()
        self.mock_os_system.assert_not_called()

    def test_single_region(self):
        """
        Tests the creation of counts for a single region.
        """
        regions = [
            ["chr10", 1000, 2000, "single_gene"]
        ]

        create_counts(self.counts_file, regions, self.vcf_compressed)

        # Verify 'open' calls
        self.mock_open.assert_has_calls([
            call(self.counts_file, "w"),
            call().__enter__(),
            call().write(""),
            call().__exit__(None, None, None),
            call(self.counts_file, "a+"),
            call().__enter__(),
            call().write("single_gene\t"),
            call().__exit__(None, None, None),
        ], any_order=False)

        # Verify os.environ calls
        region_str = "chr10:1000-2000"
        self.mock_os_environ.assert_has_calls([
            call.__setitem__('region_str', region_str),
            call.__setitem__('vcf_file', self.vcf_compressed),
            call.__setitem__('tmp_output', self.counts_file),
        ], any_order=False)
        self.assertEqual(self.mock_os_environ.__setitem__.call_count, 3)

        # Verify os.system call
        cmd = (
            'bcftools view -r "$region_str"'
            ' "$vcf_file" | grep -vc "^#" >> "$tmp_output"'
        )
        self.mock_os_system.assert_called_once_with(cmd)


class TestExtractRegions(unittest.TestCase):
    """
    Unit tests for the extract_regions function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_open = patch('builtins.open', mock_open()).start()
        self.bed_file = "/path/to/test.bed"

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_extraction_with_names(self):
        """
        Tests extraction of regions from a BED file where all regions have 
        names.
        """
        bed_content = (
            "chr1\t100\t200\tgeneA\n"
            "chr2\t300\t400\tgeneB_exon1\n"
            "chrX\t500\t600\tregion_alpha\n"
        )
        self.mock_open.return_value.__enter__.return_value = \
            bed_content.splitlines(True)

        expected_regions = [
            ("chr1", 100, 200, "geneA"),
            ("chr2", 300, 400, "geneB_exon1"),
            ("chrX", 500, 600, "region_alpha"),
        ]

        result = extract_regions(self.bed_file)

        self.mock_open.assert_called_once_with(self.bed_file)
        self.assertEqual(result, expected_regions)

    def test_successful_extraction_without_names(self):
        """
        Tests extraction of regions from a BED file where some regions 
        lack names.
        The function should generate names in the format "chr:start-end".
        """
        bed_content = (
            "chr1\t100\t200\n"
            "chr2\t300\t400\tgeneB\n"
            "chrY\t500\t600\n"
        )
        self.mock_open.return_value.__enter__.return_value = \
            bed_content.splitlines(True)

        expected_regions = [
            ("chr1", 100, 200, "chr1:100-200"),
            ("chr2", 300, 400, "geneB"),
            ("chrY", 500, 600, "chrY:500-600"),
        ]

        result = extract_regions(self.bed_file)

        self.mock_open.assert_called_once_with(self.bed_file)
        self.assertEqual(result, expected_regions)

    def test_empty_bed_file(self):
        """
        Tests extraction from an empty BED file.
        Should return an empty list.
        """
        bed_content = ""
        self.mock_open.return_value.__enter__.return_value = \
            bed_content.splitlines(True)

        expected_regions = []

        result = extract_regions(self.bed_file)

        self.mock_open.assert_called_once_with(self.bed_file)
        self.assertEqual(result, expected_regions)

    def test_bed_file_with_comments_and_empty_lines(self):
        """
        Tests extraction from a BED file containing comments and blank 
        lines.
        These lines should be ignored.
        """
        bed_content = (
            "# This is a comment line\n"
            "\n"
            "chr1\t100\t200\tgeneA\n"
            "  \t \n" # Line with only whitespace
            "chr2\t300\t400\n"
            "# Another comment\n"
            "\n"
        )
        self.mock_open.return_value.__enter__.return_value = \
            bed_content.splitlines(True)

        expected_regions = [
            ("chr1", 100, 200, "geneA"),
            ("chr2", 300, 400, "chr2:300-400"),
        ]

        result = extract_regions(self.bed_file)

        self.mock_open.assert_called_once_with(self.bed_file)
        self.assertEqual(result, expected_regions)

    def test_bed_file_with_invalid_start_end(self):
        """
        Tests extraction from a BED file with non-integer start/end values.
        This should raise a ValueError, as int() conversion will fail.
        """
        bed_content = "chr1\tinvalid_start\t200\tgeneA\n"
        self.mock_open.return_value.__enter__.return_value = \
            bed_content.splitlines(True)

        with self.assertRaises(ValueError):
            extract_regions(self.bed_file)

        self.mock_open.assert_called_once_with(self.bed_file)

    def test_bed_file_with_too_few_fields(self):
        """
        Tests extraction from a BED file with lines having fewer than 
        3 fields.
        This should raise an IndexError during tuple unpacking.
        """
        bed_content = "chr1\t100\n" # Missing 'end' field
        self.mock_open.return_value.__enter__.return_value = \
            bed_content.splitlines(True)

        with self.assertRaises(ValueError):
            extract_regions(self.bed_file)

        self.mock_open.assert_called_once_with(self.bed_file)


class TestParseAnnField(unittest.TestCase):
    """
    Unit tests for the parse_ann_field function.
    """

    def test_multiple_effects(self):
        """
        Tests parsing an ANN field with multiple effects.
        """
        ann_field = (
            "A|missense_variant|MODERATE|ENSG000001|ENST000001" +\
                "|protein_coding|1/1|HGNC:123|"
            "c.123A>G|p.Lys41Glu|123||||" +\
                "|WARNING_TRANSCRIPT_NO_START_CODON,"
            "T|stop_gained|HIGH|ENSG000002|ENST000002|" +\
                "protein_coding|1/1|HGNC:456|"
            "c.456C>T|p.Gln152*|||||"
        )
        expected_effects = ["missense_variant", "stop_gained"]
        self.assertEqual(parse_ann_field(ann_field), expected_effects)

    def test_single_effect(self):
        """
        Tests parsing an ANN field with a single effect.
        """
        ann_field = "A|frameshift_variant|HIGH|ENSG000003" +\
                "|ENST000003|protein_coding|1/1|HGNC:789" +\
                "|c.789delA|p.Met264fs|||||"
        expected_effects = ["frameshift_variant"]
        self.assertEqual(parse_ann_field(ann_field), expected_effects)

    def test_empty_ann_field(self):
        """
        Tests parsing an empty ANN field string.
        Should return an empty list.
        """
        ann_field = ""
        expected_effects = []
        self.assertEqual(parse_ann_field(ann_field), expected_effects)

    def test_ann_field_with_whitespace(self):
        """
        Tests parsing an ANN field with leading/trailing whitespace in 
        effects.
        Ensures whitespace is stripped.
        """
        ann_field = (
            "A|  synonymous_variant  |LOW|ENSG000004|ENST000004" +\
                "|protein_coding|1/1|HGNC:101|"
            "c.101G>A|p.Val34Val|||||,"
            "C| intron_variant |MODIFIER|ENSG000005|ENST000005" +\
                "|protein_coding|1/1|HGNC:112|"
            "c.112+1G>A|||||"
        )
        expected_effects = ["synonymous_variant", "intron_variant"]
        self.assertEqual(parse_ann_field(ann_field), expected_effects)

    def test_ann_field_with_malformed_entries(self):
        """
        Tests parsing an ANN field with malformed entries 
        (e.g., too few fields).
        These entries should be skipped if they don't have at least 
        two fields.
        """
        ann_field = (
            "A|missense_variant|MODERATE|ENSG000001|ENST000001" +\
                "|protein_coding|1/1|HGNC:123|"
            "c.123A>G|p.Lys41Glu|123|||||,"
            "malformed_entry_no_pipes," # This should be skipped
            "G|splice_region_variant"
        )
        expected_effects = ["missense_variant", "splice_region_variant"]
        self.assertEqual(parse_ann_field(ann_field), expected_effects)

    def test_ann_field_with_only_variant_field(self):
        """
        Tests parsing an ANN field where an annotation part only has 
        the variant field.
        Such entries should not add to effects.
        """
        # B is a malformed entry
        ann_field = "A|missense_variant|MODERATE|...,B" 
        expected_effects = ["missense_variant"]
        self.assertEqual(parse_ann_field(ann_field), expected_effects)


class TestParseGffForGenes(unittest.TestCase):
    """
    Unit tests for the parse_gff_for_genes function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_open = patch('builtins.open', mock_open()).start()
        self.gff_file_path = "/path/to/test.gff3"

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_parsing_multiple_genes(self):
        """
        Tests parsing a GFF3 file with multiple gene entries,
        including different ways gene names are specified.
        """
        gff_content = (
            "chr1\tRefSeq\tgene\t100\t500\t.\t+\t.\t" +\
                "ID=gene001;Name=GeneA;Ontology_term=GO:0003674\n"
            "chr1\tEnsembl\texon\t150\t200\t.\t+\t.\tParent=gene001\n"
            "chr2\tGnomon\tgene\t1000\t2000\t.\t-\t.\t" +\
                "ID=gene002;gene_name=GeneB;Note=hypothetical\n"
            "chrX\tCustom\tgene\t5000\t6000\t.\t+\t.\tID=gene003\n"
            "chrY\tCustom\tgene\t7000\t8000\t.\t+\t.\t" +\
                "SomeOtherAttr=value;ID=gene004;Name=GeneD\n"
        )
        self.mock_open.return_value.__enter__.return_value = \
            gff_content.splitlines(True)

        expected_regions = [
            ["chr1", 100, 500, "GeneA"],
            ["chr2", 1000, 2000, "GeneB"],
            ["chrX", 5000, 6000, "gene003"],
            ["chrY", 7000, 8000, "GeneD"],
        ]

        result = parse_gff_for_genes(self.gff_file_path)
        self.mock_open.assert_called_once_with(self.gff_file_path, 'r')
        self.assertEqual(result, expected_regions)

    def test_empty_gff_file(self):
        """
        Tests parsing an empty GFF3 file.
        Should return an empty list.
        """
        gff_content = ""
        self.mock_open.return_value.__enter__.return_value = \
            gff_content.splitlines(True)

        expected_regions = []
        result = parse_gff_for_genes(self.gff_file_path)
        self.mock_open.assert_called_once_with(self.gff_file_path, 'r')
        self.assertEqual(result, expected_regions)

    def test_gff_file_with_only_comments(self):
        """
        Tests parsing a GFF3 file containing only comment lines.
        Should return an empty list.
        """
        gff_content = (
            "##gff-version 3\n"
            "# This is a comment\n"
            "##sequence-region chr1 1 1000000\n"
        )
        self.mock_open.return_value.__enter__.return_value = \
            gff_content.splitlines(True)

        expected_regions = []
        result = parse_gff_for_genes(self.gff_file_path)
        self.mock_open.assert_called_once_with(self.gff_file_path, 'r')
        self.assertEqual(result, expected_regions)

    def test_gff_file_with_non_gene_features(self):
        """
        Tests parsing a GFF3 file containing features other than 'gene'.
        Only 'gene' features should be extracted.
        """
        gff_content = (
            "chr1\tSourceA\tgene\t100\t500\t.\t+\t.\t" +\
                "ID=gene1;Name=GeneA\n"
            "chr1\tSourceB\tmRNA\t120\t480\t.\t+\t.\tParent=gene1\n"
            "chr1\tSourceC\texon\t150\t200\t.\t+\t.\tParent=mRNA1\n"
            "chr2\tSourceD\tgene\t1000\t2000\t.\t-\t.\t" +\
                "ID=gene2;Name=GeneB\n"
            "chr2\tSourceE\tCDS\t1050\t1950\t.\t-\t0\tParent=mRNA2\n"
        )
        self.mock_open.return_value.__enter__.return_value = \
            gff_content.splitlines(True)

        expected_regions = [
            ["chr1", 100, 500, "GeneA"],
            ["chr2", 1000, 2000, "GeneB"],
        ]

        result = parse_gff_for_genes(self.gff_file_path)
        self.mock_open.assert_called_once_with(self.gff_file_path, 'r')
        self.assertEqual(result, expected_regions)

    def test_gff_file_with_malformed_lines(self):
        """
        Tests parsing a GFF3 file with lines that do not have 9 fields.
        These lines should be skipped.
        """
        gff_content = (
            "chr1\tSourceA\tgene\t100\t500\t.\t+\t.\t" +\
                "ID=gene1;Name=GeneA\n"
            "malformed_line_with_too_few_fields\n" # 1 field
            "chr2\tSourceB\tgene\t1000\t2000\t.\t-\t.\t" +\
                "ID=gene2;Name=GeneB\tExtraField\n" # 10 fields
            "chr3\tSourceC\tgene\t3000\t4000\t.\t+\t.\t" +\
                "ID=gene3;Name=GeneC\n"
        )
        self.mock_open.return_value.__enter__.return_value = \
            gff_content.splitlines(True)

        expected_regions = [
            ["chr1", 100, 500, "GeneA"],
            ["chr3", 3000, 4000, "GeneC"],
        ]

        result = parse_gff_for_genes(self.gff_file_path)
        self.mock_open.assert_called_once_with(self.gff_file_path, 'r')
        self.assertEqual(result, expected_regions)

    def test_gene_name_priority(self):
        """
        Tests the priority of gene name extraction 
        (Name > gene_name > ID > "unknown").
        """
        gff_content = (
            "chr1\tSrc\tgene\t1\t100\t.\t+\t.\t" +\
                "ID=gene_id_1;gene_name=GeneName1;" +\
                "Name=PreferredName1\n" # Name
            "chr2\tSrc\tgene\t101\t200\t.\t+\t.\t" +\
                "ID=gene_id_2;gene_name=GeneName2\n" # gene_name
            "chr3\tSrc\tgene\t201\t300\t.\t+\t.\tID=gene_id_3\n" # ID
            "chr4\tSrc\tgene\t301\t400\t.\t+\t.\t" +\
                "SomeOtherAttr=value\n" # unknown
        )
        self.mock_open.return_value.__enter__.return_value = \
            gff_content.splitlines(True)

        expected_regions = [
            ["chr1", 1, 100, "PreferredName1"],
            ["chr2", 101, 200, "GeneName2"],
            ["chr3", 201, 300, "gene_id_3"],
            ["chr4", 301, 400, "unknown"],
        ]

        result = parse_gff_for_genes(self.gff_file_path)
        self.mock_open.assert_called_once_with(self.gff_file_path, 'r')
        self.assertEqual(result, expected_regions)

    def test_gff_file_with_non_integer_coordinates(self):
        """
        Tests parsing a GFF3 file with non-integer start/end coordinates.
        This should raise a ValueError.
        """
        gff_content = "chr1\tSrc\tgene\t100\tinvalid_end" +\
                "\t.\t+\t.\tID=gene1\n"
        self.mock_open.return_value.__enter__.return_value = \
            gff_content.splitlines(True)

        with self.assertRaises(ValueError):
            parse_gff_for_genes(self.gff_file_path)

        self.mock_open.assert_called_once_with(self.gff_file_path, 'r')


class TestSamToBam(unittest.TestCase):
    """
    Unit tests for the sam_to_bam function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_mkdir = patch('os.mkdir').start()
        self.mock_system = patch('os.system').start()
        self.mock_log_code = patch('scripts.pipeline.log_code').start()

        # Define common test parameters
        self.sra_id = "test_sra_id_sam_to_bam"
        self.sam_file = "/data/input/reads.sam"
        self.output_dir = "/tmp/bam_output"
        self.log_scr = "/var/log/sam_to_bam_commands.log"

        # Pre-calculate expected file path for assertions
        self.expected_bam_file = os.path.join(
            self.output_dir,
            self.sra_id + '.bam'
        )

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_sam_to_bam_conversion(self):
        """
        Tests the full execution path of sam_to_bam when
        the output directory does not initially exist.
        Verifies all expected function calls and the return value.
        """
        # Call the function with the defined test parameters
        result_bam_file = sam_to_bam(
            sra_id=self.sra_id,
            sam_file=self.sam_file,
            output_dir=self.output_dir,
            log_scr=self.log_scr
        )

        # Assert os.mkdir is called exactly once with the correct output 
        # directory
        self.mock_mkdir.assert_called_once_with(self.output_dir)

        # Assert os.system is called with the correct samtools command
        expected_samtools_cmd = (
            f'samtools view -bS {self.sam_file}' +\
                f' -o {self.expected_bam_file}'
        )
        self.mock_system.assert_called_once_with(expected_samtools_cmd)

        # Assert log_code is called with the correct command
        self.mock_log_code.assert_called_once_with(
            expected_samtools_cmd,
            log_file=self.log_scr
        )

        # Assert the function returns the correct BAM file path
        self.assertEqual(result_bam_file, self.expected_bam_file)

    def test_sam_to_bam_directory_already_exists(self):
        """
        Tests the scenario where the output directory already exists for 
        sam_to_bam.
        Ensures that os.mkdir handles FileExistsError gracefully and
        the rest of the function proceeds as normal.
        """
        # Configure mock_mkdir to raise FileExistsError when called
        self.mock_mkdir.side_effect = FileExistsError

        # Call the function
        result_bam_file = sam_to_bam(
            sra_id=self.sra_id,
            sam_file=self.sam_file,
            output_dir=self.output_dir,
            log_scr=self.log_scr
        )

        # os.mkdir should still be attempted once, but the error is caught
        self.mock_mkdir.assert_called_once_with(self.output_dir)

        # All other calls (os.system, log_code) should still occur exactly 
        # as in the successful case.
        expected_samtools_cmd = (
            f'samtools view -bS {self.sam_file}' +\
                f' -o {self.expected_bam_file}'
        )
        self.mock_system.assert_called_once_with(expected_samtools_cmd)
        self.mock_log_code.assert_called_once_with(
            expected_samtools_cmd,
            log_file=self.log_scr
        )

        self.assertEqual(result_bam_file, self.expected_bam_file)


class TestSortIndexBam(unittest.TestCase):
    """
    Unit tests for the sort_index_bam function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_mkdir = patch('os.mkdir').start()
        self.mock_system = patch('os.system').start()

        self.mock_log_code = patch('scripts.pipeline.log_code').start()
        self.mock_log_print = patch('scripts.pipeline.log_print').start()

        # Define common test parameters for clarity and reusability
        self.sra_id = "test_sra_id_sort"
        self.bam_file = "/data/input/unsorted.bam"
        self.output_dir = "/tmp/sorted_bam_output"
        self.log_file = "/var/log/sort_pipeline.log"
        self.log_scr = "/var/log/sort_script_commands.log"

        # Pre-calculate expected file paths for assertions
        self.expected_sorted_bam = os.path.join(
            self.output_dir,
            self.sra_id + '_sorted.bam'
        )
        self.expected_bai_file = self.expected_sorted_bam + '.bai'

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_sort_index_execution(self):
        """
        Tests the full execution path of sort_index_bam when
        the output directory does not initially exist.
        Verifies all expected function calls and the return value.
        """
        # Call the function with the defined test parameters
        result_sorted_bam_file = sort_index_bam(
            sra_id=self.sra_id,
            bam_file=self.bam_file,
            output_dir=self.output_dir,
            log_file=self.log_file,
            log_scr=self.log_scr
        )

        # Assert os.mkdir is called exactly once with the correct output 
        # directory
        self.mock_mkdir.assert_called_once_with(self.output_dir)

        # Assert os.system calls in the correct order and with correct 
        # arguments
        expected_sort_cmd = f'samtools sort {self.bam_file} -o ' +\
            f'{self.expected_sorted_bam}'
        expected_index_cmd = f'samtools index {self.expected_sorted_bam}'

        # Check that os.system was called exactly 2 times
        self.assertEqual(self.mock_system.call_count, 2)
        # Verify the specific calls and their order
        self.mock_system.assert_has_calls([
            call(expected_sort_cmd),
            call(expected_index_cmd)
        ])

        # Assert log_code calls in the correct order and with correct 
        # arguments
        self.assertEqual(self.mock_log_code.call_count, 2)
        self.mock_log_code.assert_has_calls([
            call(expected_sort_cmd, log_file=self.log_scr),
            call(expected_index_cmd, log_file=self.log_scr)
        ])

        # Assert log_print calls in the correct order and with correct 
        # arguments
        self.assertEqual(self.mock_log_print.call_count, 4)
        self.mock_log_print.assert_has_calls([
            call(
                f'Sorting BAM file: {self.bam_file} ->' +\
                    f' {self.expected_sorted_bam}...',
                level='info',
                log_file=self.log_file
            ),
            call(
                f'BAM sorting complete: {self.expected_sorted_bam}',
                level='info',
                log_file=self.log_file
            ),
            call(
                f'Indexing sorted BAM file: {self.expected_sorted_bam}...',
                level='info',
                log_file=self.log_file
            ),
            call(
                f'BAM indexing complete. Index file:' +\
                    f' {self.expected_bai_file}',
                level='info',
                log_file=self.log_file
            )
        ])

        # Assert the function returns the correct sorted BAM file path
        self.assertEqual(result_sorted_bam_file, self.expected_sorted_bam)

    def test_sort_index_directory_already_exists(self):
        """
        Tests the scenario where the output directory already exists for sort_index_bam.
        Ensures that os.mkdir handles FileExistsError gracefully and
        the rest of the function proceeds as normal.
        """
        # Configure mock_mkdir to raise FileExistsError when called
        self.mock_mkdir.side_effect = FileExistsError

        # Call the function
        result_sorted_bam_file = sort_index_bam(
            sra_id=self.sra_id,
            bam_file=self.bam_file,
            output_dir=self.output_dir,
            log_file=self.log_file,
            log_scr=self.log_scr
        )

        # os.mkdir should still be attempted once, but the error is caught
        self.mock_mkdir.assert_called_once_with(self.output_dir)

        # All other calls (os.system, log_code, log_print)
        # should still occur exactly as in the successful case.
        expected_sort_cmd = f'samtools sort {self.bam_file} -o' +\
            f' {self.expected_sorted_bam}'
        expected_index_cmd = f'samtools index {self.expected_sorted_bam}'

        self.assertEqual(self.mock_system.call_count, 2)
        self.mock_system.assert_has_calls([
            call(expected_sort_cmd),
            call(expected_index_cmd)
        ])

        self.assertEqual(self.mock_log_code.call_count, 2)
        self.mock_log_code.assert_has_calls([
            call(expected_sort_cmd, log_file=self.log_scr),
            call(expected_index_cmd, log_file=self.log_scr)
        ])

        self.assertEqual(self.mock_log_print.call_count, 4)
        self.mock_log_print.assert_has_calls([
            call(
                f'Sorting BAM file: {self.bam_file} ->' +\
                    f' {self.expected_sorted_bam}...',
                level='info',
                log_file=self.log_file
            ),
            call(
                f'BAM sorting complete: {self.expected_sorted_bam}',
                level='info',
                log_file=self.log_file
            ),
            call(
                f'Indexing sorted BAM file: {self.expected_sorted_bam}...',
                level='info',
                log_file=self.log_file
            ),
            call(
                f'BAM indexing complete. Index file:' +\
                    f' {self.expected_bai_file}',
                level='info',
                log_file=self.log_file
            )
        ])

        self.assertEqual(result_sorted_bam_file, self.expected_sorted_bam)


class TestSraIdData(unittest.TestCase):
    """
    Unit tests for the sra_id_data function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_requests_get = patch('requests.get').start()
        self.mock_log_print = patch('scripts.pipeline.log_print').start()

        self.sra_id = "SRR123456"
        self.log_file = "/var/log/sra_data.log"

        # Define common expected parameters for requests.get
        self.expected_base_url = \
            "https://www.ebi.ac.uk/ena/portal/api/search"
        self.expected_fields = (
            "run_accession,fastq_ftp,sra_ftp,experiment_accession,"
            "sample_accession,study_accession,library_name," +\
                "library_strategy,"
            "library_source,library_selection,instrument_platform,"
            "instrument_model,base_count,read_count," +\
                "scientific_name,tax_id"
        )
        self.expected_params = {
            "result": "read_run",
            "query": f"run_accession={self.sra_id}",
            "fields": self.expected_fields,
            "format": "json",
            "limit": 1
        }

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_data_retrieval(self):
        """
        Tests the scenario where SRA data is successfully retrieved.
        """
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = [{
            "run_accession": self.sra_id,
            "fastq_ftp": "ftp.ebi.ac.uk/vol1/fastq/SRR123/" +\
                "SRR123456/SRR123456.fastq.gz",
            "sra_ftp": "ftp.ebi.ac.uk/vol1/sra/SRR123/" +\
                "SRR123456/SRR123456.sra",
            "experiment_accession": "ERX123456",
            "sample_accession": "ERS123456",
            "study_accession": "ERP123456",
            "library_name": "MyLib",
            "library_strategy": "WGS",
            "library_source": "GENOMIC",
            "library_selection": "RANDOM",
            "instrument_platform": "ILLUMINA",
            "instrument_model": "Illumina HiSeq 2500",
            "base_count": 1000000,
            "read_count": 500000,
            "scientific_name": "Homo sapiens",
            "tax_id": 9606
        }]
        self.mock_requests_get.return_value = mock_response

        result = sra_id_data(self.sra_id, self.log_file)

        self.mock_requests_get.assert_called_once_with(
            self.expected_base_url,
            params=self.expected_params
        )
        self.assertEqual(result, mock_response.json.return_value[0])

        # Verify log_print calls
        self.mock_log_print.assert_has_calls([
            call(
                f"Attempting to retrieve SRA data for: {self.sra_id}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f"Successfully retrieved data for {self.sra_id}.",
                level='info',
                log_file=self.log_file
            )
        ])
        self.assertEqual(self.mock_log_print.call_count, 2)

    def test_no_data_found(self):
        """
        Tests the scenario where no SRA data is found for the given ID.
        """
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = [] # Empty list indicates no data
        self.mock_requests_get.return_value = mock_response

        result = sra_id_data(self.sra_id, self.log_file)

        self.mock_requests_get.assert_called_once_with(
            self.expected_base_url,
            params=self.expected_params
        )
        self.assertIsNone(result)

        # Verify log_print calls
        self.mock_log_print.assert_has_calls([
            call(
                f"Attempting to retrieve SRA data for: {self.sra_id}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f"No SRA data found for accession ID: {self.sra_id}." +\
                    " Please check the ID.", 
                level='info',
                log_file=self.log_file
            )
        ])
        self.assertEqual(self.mock_log_print.call_count, 2)

    def test_http_error(self):
        """
        Tests handling of requests.exceptions.HTTPError.
        """
        mock_response = Mock()
        mock_response.status_code = 404
        # Configure raise_for_status to raise HTTPError
        mock_response.raise_for_status.side_effect = \
            requests.exceptions.HTTPError(
                "Not Found",
                response=mock_response
            )
        self.mock_requests_get.return_value = mock_response

        result = sra_id_data(self.sra_id, self.log_file)

        self.mock_requests_get.assert_called_once_with(
            self.expected_base_url,
            params=self.expected_params
        )
        self.assertIsNone(result)

        # Verify log_print calls
        self.mock_log_print.assert_has_calls([
            call(
                f"Attempting to retrieve SRA data for: {self.sra_id}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f"HTTP error occurred: Not Found - Status Code: 404",
                level='info',
                log_file=self.log_file
            )
        ])
        self.assertEqual(self.mock_log_print.call_count, 2)

    def test_connection_error(self):
        """
        Tests handling of requests.exceptions.ConnectionError.
        """
        self.mock_requests_get.side_effect = \
            requests.exceptions.ConnectionError("Failed to connect")

        result = sra_id_data(self.sra_id, self.log_file)

        self.mock_requests_get.assert_called_once_with(
            self.expected_base_url,
            params=self.expected_params
        )
        self.assertIsNone(result)

        # Verify log_print calls
        self.mock_log_print.assert_has_calls([
            call(
                f"Attempting to retrieve SRA data for: {self.sra_id}",
                level='info',
                log_file=self.log_file
            ),
            call(
                "Connection error occurred: Failed to connect " +\
                    "- Unable to connect to ENA API.",
                level='info',
                log_file=self.log_file
            )
        ])
        self.assertEqual(self.mock_log_print.call_count, 2)

    def test_timeout_error(self):
        """
        Tests handling of requests.exceptions.Timeout.
        """
        self.mock_requests_get.side_effect = \
            requests.exceptions.Timeout("Request timed out")

        result = sra_id_data(self.sra_id, self.log_file)

        self.mock_requests_get.assert_called_once_with(
            self.expected_base_url,
            params=self.expected_params
        )
        self.assertIsNone(result)

        # Verify log_print calls
        self.mock_log_print.assert_has_calls([
            call(
                f"Attempting to retrieve SRA data for: {self.sra_id}",
                level='info',
                log_file=self.log_file
            ),
            call(
                "Timeout error occurred: Request timed out " +\
                    "- Request to ENA API timed out.",
                level='info',
                log_file=self.log_file
            )
        ])
        self.assertEqual(self.mock_log_print.call_count, 2)

    def test_request_exception(self):
        """
        Tests handling of a general requests.exceptions.RequestException.
        """
        self.mock_requests_get.side_effect = \
            requests.exceptions.RequestException("Generic request error")

        result = sra_id_data(self.sra_id, self.log_file)

        self.mock_requests_get.assert_called_once_with(
            self.expected_base_url,
            params=self.expected_params
        )
        self.assertIsNone(result)

        # Verify log_print calls
        self.mock_log_print.assert_has_calls([
            call(
                f"Attempting to retrieve SRA data for: {self.sra_id}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f"An unexpected error occurred during the request: " +\
                    "Generic request error",
                level='info',
                log_file=self.log_file
            )
        ])
        self.assertEqual(self.mock_log_print.call_count, 2)

    def test_json_decode_error(self):
        """
        Tests handling of json.JSONDecodeError when response content 
        is invalid JSON.
        """
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.side_effect = json.JSONDecodeError(
            "Invalid JSON",
            doc="malformed",
            pos=0
        )
        mock_response.text = "This is not JSON" # Raw response text
        self.mock_requests_get.return_value = mock_response

        result = sra_id_data(self.sra_id, self.log_file)

        self.mock_requests_get.assert_called_once_with(
            self.expected_base_url,
            params=self.expected_params
        )
        self.assertIsNone(result)

        # Verify log_print calls
        self.mock_log_print.assert_has_calls([
            call(
                f"Attempting to retrieve SRA data for: {self.sra_id}",
                level='info',
                log_file=self.log_file
            ),
            call(
                "Error decoding JSON response: Invalid JSON: " +\
                    "line 1 column 1 (char 0). Response content: " +\
                    "This is not JSON",
                level='info',
                log_file=self.log_file
            )
        ])
        self.assertEqual(self.mock_log_print.call_count, 2)

    def test_generic_exception(self):
        """
        Tests handling of any other unexpected Exception.
        """
        self.mock_requests_get.side_effect = Exception(
            "Something went wrong unexpectedly"
        )

        result = sra_id_data(self.sra_id, self.log_file)

        self.mock_requests_get.assert_called_once_with(
            self.expected_base_url,
            params=self.expected_params
        )
        self.assertIsNone(result)

        # Verify log_print calls 
        # (level should be 'warn' for generic Exception)
        self.mock_log_print.assert_has_calls([
            call(
                f"Attempting to retrieve SRA data for: {self.sra_id}",
                level='info',
                log_file=self.log_file
            ),
            call(
                "An unexpected error occurred: " +\
                    "Something went wrong unexpectedly",
                level='warn',
                log_file=self.log_file
            )
        ])
        self.assertEqual(self.mock_log_print.call_count, 2)


class TestSraToFastq(unittest.TestCase):
    """
    Unit tests for the sra_to_fastq function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_mkdir = patch('os.mkdir').start()
        self.mock_system = patch('os.system').start()
        self.mock_exists = patch('os.path.exists').start()
        self.mock_remove = patch('os.remove').start()
        self.mock_log_print = patch('scripts.pipeline.log_print').start()
        self.mock_log_code = patch('scripts.pipeline.log_code').start()
        self.mock_time_sleep = patch('time.sleep').start()

        self.sra_id = "test_sra_id_fastq"
        self.output_dir = "/tmp/fastq_output"
        self.log_file = "/var/log/fastq_pipeline.log"
        self.log_scr = "/var/log/fastq_script.log"

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_download_single_fastq(self):
        l_fastq = [os.path.join(
            self.output_dir,
            f"{self.sra_id}.fastq.gz"
        )]
        # Simulate success
        self.mock_system.return_value = 0
        # Files don't exist before download
        self.mock_exists.return_value = False

        sra_to_fastq(
            self.sra_id,
            l_fastq,
            self.output_dir,
            self.log_file,
            self.log_scr
        )

        self.mock_mkdir.assert_called_once_with(self.output_dir)
        expected_cmd = f'fastq-dump --gzip -O {self.output_dir}' +\
            f' {self.sra_id}'
        self.mock_log_code.assert_called_once_with(
            expected_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_called_once_with(expected_cmd)
        # No errors, no warnings
        self.mock_log_print.assert_not_called()
        self.mock_time_sleep.assert_not_called()
        self.mock_remove.assert_not_called()

    def test_successful_download_paired_fastq(self):
        l_fastq = [
            os.path.join(self.output_dir, f"{self.sra_id}_1.fastq.gz"),
            os.path.join(self.output_dir, f"{self.sra_id}_2.fastq.gz")
        ]
        self.mock_system.return_value = 0 # Simulate success
        self.mock_exists.return_value = False

        sra_to_fastq(
            self.sra_id,
            l_fastq,
            self.output_dir,
            self.log_file,
            self.log_scr
        )

        self.mock_mkdir.assert_called_once_with(self.output_dir)
        expected_cmd = f'fastq-dump --gzip --split-files ' +\
            f'-O {self.output_dir} {self.sra_id}'
        self.mock_log_code.assert_called_once_with(
            expected_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_called_once_with(expected_cmd)
        self.mock_log_print.assert_not_called()
        self.mock_time_sleep.assert_not_called()
        self.mock_remove.assert_not_called()

    def test_sra_to_fastq_directory_already_exists(self):
        self.mock_mkdir.side_effect = FileExistsError
        l_fastq = [os.path.join(
            self.output_dir,
            f"{self.sra_id}.fastq.gz"
        )]
        self.mock_system.return_value = 0
        self.mock_exists.return_value = False

        sra_to_fastq(
            self.sra_id,
            l_fastq,
            self.output_dir,
            self.log_file,
            self.log_scr
        )
        self.mock_mkdir.assert_called_once_with(self.output_dir)
        expected_cmd = f'fastq-dump --gzip -O {self.output_dir}' +\
            f' {self.sra_id}'
        self.mock_system.assert_called_once_with(expected_cmd)

    def test_fastq_dump_failure_then_success(self):
        l_fastq = [os.path.join(
            self.output_dir,
            f"{self.sra_id}.fastq.gz"
        )]
        # Simulate failure on first attempt 
        # (os.system returns non-zero), then success
        self.mock_system.side_effect = [Exception('fail'), 0]
        # Simulate files not existing on first check, then existing 
        # for cleanup, then not existing for next attempt
        self.mock_exists.side_effect = [True, False]
        self.mock_remove.return_value = None # os.remove succeeds

        sra_to_fastq(
            self.sra_id,
            l_fastq,
            self.output_dir,
            self.log_file,
            self.log_scr,
            max_retry=2
        )

        expected_cmd = f'fastq-dump --gzip -O {self.output_dir}' +\
            f' {self.sra_id}'
        self.assertEqual(self.mock_system.call_count, 2)
        self.mock_system.assert_has_calls(
            [call(expected_cmd), call(expected_cmd)]
        )

        self.mock_log_print.assert_called_once_with(
            mock_ANY, # Message will contain "Attempt 1 failed"
            level='warn',
            log_file=self.log_file
        )
        self.assertIn(
            "Attempt 1 failed",
            self.mock_log_print.call_args[0][0]
        )
        self.mock_time_sleep.assert_called_once_with(60)
        self.mock_remove.assert_called_once_with(l_fastq[0])

    def test_fastq_dump_persistent_failure(self):
        l_fastq = [os.path.join(
            self.output_dir,
            f"{self.sra_id}.fastq.gz"
        )]
        # Define max_retry
        max_retry = 2
        # Simulate consistent failure for all retries
        self.mock_system.side_effect = [Exception('fail')] * max_retry
        # Simulate files existing for cleanup after each failed attempt
        self.mock_exists.return_value = True
        self.mock_remove.return_value = None

        sra_to_fastq(
            self.sra_id,
            l_fastq,
            self.output_dir,
            self.log_file,
            self.log_scr,
            max_retry=max_retry
        )

        expected_cmd = f'fastq-dump --gzip -O {self.output_dir}' +\
            f' {self.sra_id}'
        self.assertEqual(self.mock_system.call_count, max_retry)
        self.mock_system.assert_has_calls(
            [call(expected_cmd), call(expected_cmd)]
        )

        self.assertEqual(self.mock_log_print.call_count, max_retry+1)
        self.mock_log_print.assert_any_call(
            mock_ANY,
            level='warn',
            log_file=self.log_file
        ) # Attempt 1 failed
        self.mock_log_print.assert_any_call(
            mock_ANY,
            level='warn',
            log_file=self.log_file
        ) 
        # Attempt 2 to max_retry failed
        for re in range(1, max_retry):
            self.mock_log_print.assert_any_call(
                f"Failed to download FASTQ files after {re+1} attempts.",
                level='error',
                log_file=self.log_file
            )
        # Only sleeps after first failed attempt
        self.assertEqual(self.mock_time_sleep.call_count, max_retry-1)
        # remove called after each failed attempt
        self.assertEqual(self.mock_remove.call_count, max_retry)

    def test_fastq_dump_exception_during_run(self):
        l_fastq = [os.path.join(
            self.output_dir,
            f"{self.sra_id}.fastq.gz"
        )]
        # Simulate an exception (e.g., permission error) during 
        # os.system call
        self.mock_system.side_effect = Exception("Mock system error")
        self.mock_exists.return_value = False # No files to clean up

        sra_to_fastq(
            self.sra_id,
            l_fastq,
            self.output_dir,
            self.log_file,
            self.log_scr,
            max_retry=1
        )

        expected_cmd = f'fastq-dump --gzip -O {self.output_dir}' +\
            f' {self.sra_id}'
        self.mock_system.assert_called_once_with(expected_cmd)

        self.assertEqual(self.mock_log_print.call_count, 2)
        self.mock_log_print.assert_any_call(
            f"Attempt 1 failed: Mock system error",
            level='warn',
            log_file=self.log_file
        )
        self.mock_log_print.assert_any_call(
            f"Failed to download FASTQ files after 1 attempts.",
            level='error',
            log_file=self.log_file
        )
        self.mock_time_sleep.assert_not_called()
        self.mock_remove.assert_not_called()


class TestVariantAnalysisSnpEff(unittest.TestCase):
    """
    Unit tests for the variant_analysis_snpeff function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_mkdir = patch('os.mkdir').start()
        self.mock_system = patch('os.system').start()
        self.mock_log_print = patch('scripts.pipeline.log_print').start()
        self.mock_log_code = patch('scripts.pipeline.log_code').start()

        self.sra_id = "test_sra_id_snpeff"
        self.vcf_file = "/data/input.vcf"
        self.genome_name = "GRCh38"
        self.snpeff_dir = "/opt/snpEff"
        self.output_dir = "/tmp/snpeff_output"
        self.log_file = "/var/log/snpeff_pipeline.log"
        self.log_scr = "/var/log/snpeff_script.log"
        self.expected_snpeff_vcf = os.path.join(
            self.output_dir,
            self.sra_id + '_snpeff.vcf'
        )
        self.expected_snpeff_vcf_gz = self.expected_snpeff_vcf + '.gz'
        self.expected_snpeff_vcf_tbi = self.expected_snpeff_vcf_gz + '.tbi'

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_snpeff_analysis(self):
        # Simulate success for all commands
        self.mock_system.return_value = 0

        result = variant_analysis_snpeff(
            self.sra_id,
            self.vcf_file,
            self.genome_name,
            self.snpeff_dir,
            self.output_dir,
            self.log_file,
            self.log_scr
        )

        self.mock_mkdir.assert_called_once_with(self.output_dir)

        # Assert snpEff command and logs
        expected_snpeff_cmd = (
            f'java -Xmx8g -jar {self.snpeff_dir}/snpEff/snpEff.jar '
            f'{self.genome_name} {self.vcf_file} '
            f'> {self.expected_snpeff_vcf}'
        )
        self.mock_log_print.assert_any_call(
            f"Analyzing variants from {self.vcf_file} with snpEff...",
            level='info',
            log_file=self.log_file
        )
        self.mock_log_code.assert_any_call(
            expected_snpeff_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_any_call(expected_snpeff_cmd)

        # Assert tail command
        expected_tail_cmd = f'tail {self.expected_snpeff_vcf}'
        self.mock_log_code.assert_any_call(
            expected_tail_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_any_call(expected_tail_cmd)

        # Assert bgzip command and logs
        expected_bgzip_cmd = f'bgzip -c {self.expected_snpeff_vcf} ' +\
            f'> {self.expected_snpeff_vcf_gz}'
        self.mock_log_print.assert_any_call(
            f"Compressing {self.expected_snpeff_vcf} with bgzip...",
            level='info',
            log_file=self.log_file
        )
        self.mock_log_code.assert_any_call(
            expected_bgzip_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_any_call(expected_bgzip_cmd)
        self.mock_log_print.assert_any_call(
            f"Compressed VCF: {self.expected_snpeff_vcf_gz}",
            level='info',
            log_file=self.log_file
        )

        # Assert tabix command and logs
        expected_tabix_cmd = f'tabix -p vcf {self.expected_snpeff_vcf_gz}'
        self.mock_log_print.assert_any_call(
            f"Indexing {self.expected_snpeff_vcf_gz} with tabix...",
            level='info',
            log_file=self.log_file
        )
        self.mock_log_code.assert_any_call(
            expected_tabix_cmd,
            log_file=self.log_scr
        )
        self.mock_system.assert_any_call(expected_tabix_cmd)
        self.mock_log_print.assert_any_call(
            f"VCF index: {self.expected_snpeff_vcf_tbi}",
            level='info',
            log_file=self.log_file
        )

        self.assertEqual(result, self.expected_snpeff_vcf_gz)
        self.assertEqual(self.mock_system.call_count, 4)
        self.assertEqual(self.mock_log_print.call_count, 5)
        self.assertEqual(self.mock_log_code.call_count, 4)

    def test_snpeff_analysis_directory_already_exists(self):
        self.mock_mkdir.side_effect = FileExistsError
        self.mock_system.return_value = 0

        variant_analysis_snpeff(
            self.sra_id,
            self.vcf_file,
            self.genome_name,
            self.snpeff_dir,
            self.output_dir,
            self.log_file,
            self.log_scr
        )
        self.mock_mkdir.assert_called_once_with(self.output_dir)
        # All commands should still run
        self.assertEqual(self.mock_system.call_count, 4)

    def test_tabix_command_failure(self):
        # snpEff, tail, bgzip success; tabix fails
        self.mock_system.side_effect = [0, 0, 0, 1]

        result = variant_analysis_snpeff(
            self.sra_id,
            self.vcf_file,
            self.genome_name,
            self.snpeff_dir,
            self.output_dir,
            self.log_file,
            self.log_scr
        )

        expected_tabix_cmd = f'tabix -p vcf {self.expected_snpeff_vcf_gz}'
        self.mock_system.assert_has_calls([
            call(mock_ANY), # snpEff
            call(mock_ANY), # tail
            call(mock_ANY), # bgzip
            call(expected_tabix_cmd)
        ])
        self.assertEqual(self.mock_system.call_count, 4)
        self.assertEqual(result, self.expected_snpeff_vcf_gz)
        self.assertEqual(self.mock_log_print.call_count, 5)
        self.assertEqual(self.mock_log_code.call_count, 4)


class TestVariantCallMpileup(unittest.TestCase):
    """
    Unit tests for the variant_call_mpileup function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_mkdir = patch('os.mkdir').start()
        self.mock_system = patch('os.system').start()
        self.mock_remove = patch('os.remove').start()

        self.mock_log_code = patch('scripts.pipeline.log_code').start()
        self.mock_log_print = patch('scripts.pipeline.log_print').start()

        # Define common test parameters for clarity and reusability
        self.sra_id = "test_sra_id"
        self.sorted_bam = "/data/input/aligned.bam"
        self.ref_genome = "/data/references/genome.fasta"
        self.output_dir = "/tmp/variant_output"
        self.log_file = "/var/log/pipeline.log"
        self.log_scr = "/var/log/script_commands.log"

        # Pre-calculate expected file paths for assertions
        self.expected_bcf_file = os.path.join(
            self.output_dir,
            self.sra_id + '.bcf'
        )
        self.expected_vcf_file = os.path.join(
            self.output_dir,
            self.sra_id + '.vcf'
        )

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_variant_call_mpileup(self):
        """
        Tests the full execution of variant_call_mpileup when
        the output directory does not initially exist.
        Verifies all expected function calls and the return value.
        """
        # Call the function with the defined test parameters
        result_vcf_file = variant_call_mpileup(
            sra_id=self.sra_id,
            sorted_bam=self.sorted_bam,
            ref_genome=self.ref_genome,
            output_dir=self.output_dir,
            log_file=self.log_file,
            log_scr=self.log_scr
        )

        # Assert os.mkdir is called exactly once with the correct output 
        # directory
        self.mock_mkdir.assert_called_once_with(self.output_dir)

        # Assert os.system calls in the correct order and with correct 
        # arguments
        expected_mpileup_cmd = (
            f'bcftools mpileup -f {self.ref_genome}' +\
                f' {self.sorted_bam} > {self.expected_bcf_file}'
        )
        expected_tail_bcf_cmd = f'tail {self.expected_bcf_file}'
        expected_call_cmd = (
            f'bcftools call -mv -o {self.expected_vcf_file}' +\
                f' {self.expected_bcf_file}'
        )
        expected_tail_vcf_cmd = f'tail {self.expected_vcf_file}'

        # Check that os.system was called exactly 4 times
        self.assertEqual(self.mock_system.call_count, 4)
        # Verify the specific calls and their order
        self.mock_system.assert_has_calls([
            call(expected_mpileup_cmd),
            call(expected_tail_bcf_cmd),
            call(expected_call_cmd),
            call(expected_tail_vcf_cmd)
        ])

        # Assert os.remove is called exactly once with the intermediate 
        # BCF file
        self.mock_remove.assert_called_once_with(self.expected_bcf_file)

        # Assert log_code calls in the correct order and with correct 
        # arguments
        self.assertEqual(self.mock_log_code.call_count, 4)
        self.mock_log_code.assert_has_calls([
            call(expected_mpileup_cmd, log_file=self.log_scr),
            call(expected_tail_bcf_cmd, log_file=self.log_scr),
            call(expected_call_cmd, log_file=self.log_scr),
            call(expected_tail_vcf_cmd, log_file=self.log_scr)
        ])

        # Assert log_print calls in the correct order and with correct 
        # arguments
        self.assertEqual(self.mock_log_print.call_count, 2)
        self.mock_log_print.assert_has_calls([
            call(
                f"Pileup and BCF file generated: {self.expected_bcf_file}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f'Variant calling complete. Output VCF:' +\
                f' {self.expected_vcf_file}',
                level='info',
                log_file=self.log_file
            )
        ])

        # Assert the function returns the correct VCF file path
        self.assertEqual(result_vcf_file, self.expected_vcf_file)

    def test_variant_call_mpileup_directory_exists(self):
        """
        Tests the scenario where the output directory already exists.
        Ensures that os.mkdir handles FileExistsError gracefully and
        the rest of the function proceeds as normal.
        """
        # Configure mock_mkdir to raise FileExistsError when called
        self.mock_mkdir.side_effect = FileExistsError

        # Call the function
        result_vcf_file = variant_call_mpileup(
            sra_id=self.sra_id,
            sorted_bam=self.sorted_bam,
            ref_genome=self.ref_genome,
            output_dir=self.output_dir,
            log_file=self.log_file,
            log_scr=self.log_scr
        )

        # os.mkdir should still be attempted once, but the error is caught
        self.mock_mkdir.assert_called_once_with(self.output_dir)

        # All other calls (os.system, os.remove, log_code, log_print)
        # should still occur exactly as in the successful case,
        # as the FileExistsError is handled and does not stop execution.
        expected_mpileup_cmd = (
            f'bcftools mpileup -f {self.ref_genome}' +\
                f' {self.sorted_bam} > {self.expected_bcf_file}'
        )
        expected_tail_bcf_cmd = f'tail {self.expected_bcf_file}'
        expected_call_cmd = (
            f'bcftools call -mv -o {self.expected_vcf_file}' +\
                f' {self.expected_bcf_file}'
        )
        expected_tail_vcf_cmd = f'tail {self.expected_vcf_file}'

        self.assertEqual(self.mock_system.call_count, 4)
        self.mock_system.assert_has_calls([
            call(expected_mpileup_cmd),
            call(expected_tail_bcf_cmd),
            call(expected_call_cmd),
            call(expected_tail_vcf_cmd)
        ])
        self.mock_remove.assert_called_once_with(self.expected_bcf_file)

        self.assertEqual(self.mock_log_code.call_count, 4)
        self.mock_log_code.assert_has_calls([
            call(expected_mpileup_cmd, log_file=self.log_scr),
            call(expected_tail_bcf_cmd, log_file=self.log_scr),
            call(expected_call_cmd, log_file=self.log_scr),
            call(expected_tail_vcf_cmd, log_file=self.log_scr)
        ])

        self.assertEqual(self.mock_log_print.call_count, 2)
        self.mock_log_print.assert_has_calls([
            call(
                f"Pileup and BCF file generated: {self.expected_bcf_file}",
                level='info',
                log_file=self.log_file
            ),
            call(
                f'Variant calling complete. Output VCF:' +\
                    f' {self.expected_vcf_file}',
                level='info',
                log_file=self.log_file
            )
        ])

        self.assertEqual(result_vcf_file, self.expected_vcf_file)


class TestVariantsPerBinOs(unittest.TestCase):
    """
    Unit tests for the variants_per_bin_os function.
    """

    def setUp(self):
        """
        Set up mocks before each test method runs.
        """
        self.mock_tempfile_tmpdir = patch(
            'tempfile.TemporaryDirectory'
        ).start()
        self.mock_os_system = patch('os.system').start()
        self.mock_gzip_open = patch('gzip.open', mock_open()).start()
        self.mock_builtins_open = patch(
            'builtins.open',
            mock_open()
        ).start()
        self.mock_pd_read_csv = patch('pandas.read_csv').start()
        self.mock_log_code = patch('scripts.pipeline.log_code').start()
        self.mock_log_print = patch('scripts.pipeline.log_print').start()

        # Define common test parameters
        self.vcf_file = "/path/to/variants.vcf.gz"
        self.genome_sizes = "/path/to/genome.sizes"
        self.bin_size = 10000
        self.log_file = "/var/log/pipeline.log"
        self.log_scr = "/var/log/script_commands.log"

        # Simulate a temporary directory path
        self.tmpdir_path = "/mock/tmp/dir"
        self.mock_tempfile_tmpdir.return_value.__enter__.return_value = \
            self.tmpdir_path

        # Expected intermediate file paths
        self.bins_bed = os.path.join(self.tmpdir_path, "genome_bins.bed")
        self.variants_bed = os.path.join(self.tmpdir_path, "variants.bed")
        self.intersected = os.path.join(
            self.tmpdir_path,
            "intersected.bed"
        )

    def tearDown(self):
        """
        Clean up mocks after each test method runs.
        """
        patch.stopall()

    def test_successful_variants_per_bin_counting(self):
        """
        Tests the full successful execution path of variants_per_bin_os.
        """
        # Mock bedtools makewindows to succeed
        self.mock_os_system.side_effect = [0, 0]

        # Mock VCF content for gzip.open
        mock_vcf_content = [
            "##header_line\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t" +\
                "INFO\tFORMAT\tSAMPLE\n",
            "chr1\t101\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            "chr1\t20001\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\n",
            "chr2\t500\t.\tAT\tA\t.\tPASS\t.\tGT\t0/0\n"
        ]
        self.mock_gzip_open.return_value.__enter__.return_value = \
            mock_vcf_content

        # Mock pandas.read_csv to return a DataFrame
        mock_df_data = {
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [0, 10000, 0],
            "end": [10000, 20000, 10000],
            "variant_count": [1, 1, 1]
        }
        self.mock_pd_read_csv.return_value = pd.DataFrame(mock_df_data)

        expected_output = {
            "variants_in_chr1:0-10000": 1,
            "variants_in_chr1:10000-20000": 1,
            "variants_in_chr2:0-10000": 1
        }

        result = variants_per_bin_os(
            self.vcf_file,
            self.genome_sizes,
            self.bin_size,
            self.log_file,
            self.log_scr
        )

        # Assert bedtools makewindows command
        expected_makewindows_cmd = "bedtools makewindows " +\
            f"-g {self.genome_sizes} -w {self.bin_size} > {self.bins_bed}"
        self.mock_log_code.assert_any_call(
            expected_makewindows_cmd,
            log_file=self.log_scr
        )
        self.mock_os_system.assert_any_call(expected_makewindows_cmd)

        # Assert VCF to BED conversion file operations
        self.mock_gzip_open.assert_called_once_with(self.vcf_file, 'rt')
        self.mock_builtins_open.assert_any_call(self.variants_bed, "w")
        # Check content written to variants_bed
        self.mock_builtins_open.return_value.__enter__.return_value.write.assert_has_calls([
            call("chr1\t100\t101\n"),
            call("chr1\t20000\t20001\n"),
            call("chr2\t499\t501\n")
        ])

        # Assert bedtools intersect command
        expected_intersect_cmd = "bedtools intersect " +\
            f"-a {self.bins_bed} -b {self.variants_bed} " +\
            f"-c > {self.intersected}"
        self.mock_log_code.assert_any_call(
            expected_intersect_cmd,
            log_file=self.log_scr
        )
        self.mock_os_system.assert_any_call(expected_intersect_cmd)

        # Assert pandas.read_csv call
        self.mock_pd_read_csv.assert_called_once_with(
            self.intersected,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "variant_count"]
        )

        # Assert final output dictionary
        self.assertEqual(result, expected_output)

        # Assert log_print calls
        self.mock_log_print.assert_not_called()

    def test_makewindows_command_failure(self):
        """
        Tests the scenario where bedtools makewindows command fails.
        """
        self.mock_os_system.return_value = 1 # Simulate failure

        result = variants_per_bin_os(
            self.vcf_file,
            self.genome_sizes,
            self.bin_size,
            self.log_file,
            self.log_scr
        )

        expected_makewindows_cmd = "bedtools makewindows " +\
            f"-g {self.genome_sizes} -w {self.bin_size} > {self.bins_bed}"
        self.mock_os_system.assert_called_once_with(
            expected_makewindows_cmd
        )
        self.assertEqual(result, {})
        self.mock_log_print.assert_called_once_with(
            "Failed to generate genome bins using bedtools.",
            level='warn',
            log_file=self.log_file
        )
        # Should not proceed to VCF conversion
        self.mock_gzip_open.assert_not_called()
        self.mock_pd_read_csv.assert_not_called()

    def test_vcf_to_bed_conversion_exception(self):
        """
        Tests handling of an exception during VCF to BED conversion 
        (e.g., gzip.open error).
        """
        self.mock_os_system.return_value = 0 # makewindows success
        self.mock_gzip_open.side_effect = Exception(
            "Mock gzip open error"
        )

        result = variants_per_bin_os(
            self.vcf_file,
            self.genome_sizes,
            self.bin_size,
            self.log_file,
            self.log_scr
        )

        self.assertEqual(result, {})
        self.mock_log_print.assert_called_once_with(
            "Failed to convert VCF to BED: Mock gzip open error",
            level='warn',
            log_file=self.log_file
        )
        # Only makewindows should have run
        self.mock_os_system.assert_called_once()
        self.mock_pd_read_csv.assert_not_called()

    def test_vcf_to_bed_conversion_malformed_vcf_data(self):
        """
        Tests handling of malformed VCF data that causes errors during 
        parsing.
        """
        self.mock_os_system.return_value = 0 # makewindows success
        mock_vcf_content = [
            "##header\n",
            "chr1\tPOS\tID\tREF\tALT\n", # Malformed line (POS not int)
        ]
        self.mock_gzip_open.return_value.__enter__.return_value = \
            mock_vcf_content

        result = variants_per_bin_os(
            self.vcf_file,
            self.genome_sizes,
            self.bin_size,
            self.log_file,
            self.log_scr
        )

        self.assertEqual(result, {})
        # Check that the log_print for ValueError from int() conversion 
        # is called
        self.mock_log_print.assert_called_once_with(
            mock_ANY, # Message will contain ValueError details
            level='warn',
            log_file=self.log_file
        )
        self.assertIn(
            "Failed to convert VCF to BED: invalid literal for int()",
            self.mock_log_print.call_args[0][0]
        )
        # Only makewindows should have run
        self.mock_os_system.assert_called_once()
        self.mock_pd_read_csv.assert_not_called()


    def test_intersect_command_failure(self):
        """
        Tests the scenario where bedtools intersect command fails.
        """
        self.mock_os_system.side_effect = [0, 1]

        # Mock VCF content for gzip.open
        mock_vcf_content = ["chr1\t101\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"]
        self.mock_gzip_open.return_value.__enter__.return_value = \
            mock_vcf_content

        result = variants_per_bin_os(
            self.vcf_file,
            self.genome_sizes,
            self.bin_size,
            self.log_file,
            self.log_scr
        )

        expected_makewindows_cmd = f"bedtools makewindows " +\
            f"-g {self.genome_sizes} -w {self.bin_size} > {self.bins_bed}"
        expected_intersect_cmd = f"bedtools intersect " +\
            f"-a {self.bins_bed} -b {self.variants_bed} " +\
            f"-c > {self.intersected}"

        self.mock_os_system.assert_has_calls([
            call(expected_makewindows_cmd),
            call(expected_intersect_cmd)
        ])
        self.assertEqual(result, {})
        self.mock_log_print.assert_called_once_with(
            "Failed to intersect variants with genome bins.",
            level='warn',
            log_file=self.log_file
        )
        # Should not proceed to read_csv
        self.mock_pd_read_csv.assert_not_called()

    def test_read_csv_exception(self):
        """
        Tests handling of an exception during pandas.read_csv.
        """
        self.mock_os_system.side_effect = [0, 0]

        # Mock VCF content for gzip.open
        mock_vcf_content = ["chr1\t101\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"]
        self.mock_gzip_open.return_value.__enter__.return_value = \
            mock_vcf_content

        self.mock_pd_read_csv.side_effect = Exception(
            "Mock pandas read_csv error"
        )

        result = variants_per_bin_os(
            self.vcf_file,
            self.genome_sizes,
            self.bin_size,
            self.log_file,
            self.log_scr
        )

        self.assertEqual(result, {})
        self.mock_log_print.assert_called_once_with(
            "Failed to read intersection result: " +\
                "Mock pandas read_csv error",
            level='warn',
            log_file=self.log_file
        )
        # Both bedtools commands should have run
        self.mock_os_system.call_count == 2

    def test_empty_intersected_data(self):
        """
        Tests the scenario where the intersected BED file is empty,
        resulting in an empty DataFrame.
        """
        self.mock_os_system.side_effect = [0, 0]

        # Mock VCF content for gzip.open
        mock_vcf_content = ["chr1\t101\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"]
        self.mock_gzip_open.return_value.__enter__.return_value = \
            mock_vcf_content

        # Mock pandas.read_csv to return an empty DataFrame
        self.mock_pd_read_csv.return_value = pd.DataFrame(
            columns=["chrom", "start", "end", "variant_count"]
        )

        result = variants_per_bin_os(
            self.vcf_file,
            self.genome_sizes,
            self.bin_size,
            self.log_file,
            self.log_scr
        )

        self.assertEqual(result, {})
        self.mock_log_print.assert_called_once_with(
            "No intersected variants found.",
            level='warn',
            log_file=self.log_file
        )
        # Ensure bedtools commands and VCF conversion still ran
        self.mock_os_system.call_count == 2
        self.mock_gzip_open.assert_called_once()
        self.mock_pd_read_csv.assert_called_once()


if __name__=='__main__':
    unittest.main()
