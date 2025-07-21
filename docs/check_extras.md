# check_extras.R

## Purpose and Role
This file contains the **parameter validation and configuration system** for the polyester RNA-seq simulation pipeline. It provides the `.check_extras()` function that validates, sanitizes, and sets intelligent defaults for all optional parameters passed via the `...` argument in simulation functions, ensuring robust parameter handling across different simulation scenarios.

## Pipeline Flow
1. **Parameter Reception**: Receives `extras` list from `...` arguments
2. **Validation**: Validates each parameter against allowed values and constraints
3. **Default Setting**: Assigns sensible defaults for missing parameters
4. **Compatibility Checking**: Ensures parameter combinations are compatible
5. **Sanitization**: Cleans and formats parameter values for downstream use
6. **Return**: Returns validated and populated `extras` list

## Key Variables

### Input Parameters (in `extras` list)
- `distr`: Fragment length distribution type ('normal', 'empirical', 'custom')
- `fraglen`: Mean fragment length (default: 250)
- `fragsd`: Fragment length standard deviation (default: 25)
- `readlen`: Read length (default: 100)
- `bias`: Positional bias model ('none', 'rnaf', 'cdnaf')
- `error_model`: Error model ('uniform', 'illumina4', 'illumina5', 'custom')
- `error_rate`: Overall error rate for uniform model (default: 0.005)
- `lib_sizes`: Library size factors for each sample (default: 1)
- `frag_GC_bias`: GC bias matrix or 'none' (default: 'none')
- `strand_specific`: Strand-specific simulation flag (default: FALSE)
- `shuffle`: Read shuffling flag (default: FALSE)
- `fastq`: FASTQ output flag (default: FALSE)
- `verbose`: Progress reporting flag (default: FALSE)
- `seq_depth`: Target sequencing depth per sample
- `pcr_rate`: PCR duplication rate (default: NULL)
- `pcr_lambda`: PCR duplication lambda (default: 1)
- `adapter_contamination`: Adapter contamination flag (default: FALSE)
- `adapter_sequence`: Adapter sequence for contamination (default: 'CTGTCTCTTATACACATCT')

### Validation Functions
- `.check_error_model()`: Validates error model parameters
- `.check_fold_changes()`: Validates fold change matrix dimensions
- `.write_info()`: Writes simulation information files

## Complex Structures

### Parameter Validation Rules
- **Distribution Validation**: Ensures `distr` is one of allowed values
- **Numeric Constraints**: Validates ranges for error rates, fragment lengths, etc.
- **Matrix Dimensions**: Checks GC bias matrix dimensions (101 × samples)
- **Compatibility Checks**: Ensures custom error models have required files
- **Length Matching**: Validates parameter vectors match sample count

### Default Assignment Logic
- **Smart Defaults**: Uses sensible defaults based on typical RNA-seq parameters
- **Vector Expansion**: Expands single values to match sample count
- **Conditional Defaults**: Sets defaults based on other parameter values
- **Type Coercion**: Ensures consistent data types across parameters

## Key Dependencies
- **File System**: Checks existence of custom error model files
- **Data Types**: Uses `data.table` for efficient matrix operations
- **Validation**: Provides comprehensive parameter checking across all simulation functions

## Usage Notes
This validation system acts as a central hub for parameter handling, ensuring that all simulation functions (simulate_experiment, simulate_experiment_countmat, simulate_experiment_empirical) use consistent parameter validation and default assignment logic.