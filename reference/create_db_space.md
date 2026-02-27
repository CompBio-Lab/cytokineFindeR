# Create an LRI database space

Create an LRI database space

## Usage

``` r
create_db_space(
  ...,
  filePath = "data/ligand_receptor_db.rda",
  saveToFile = TRUE
)
```

## Arguments

- ...:

  Named database objects to save (passed as name = object pairs).

- filePath:

  Character string specifying the file path for the saved `.rda` file.
  Defaults to `"data/ligand_receptor_db.rda"`.

- saveToFile:

  Logical. If `TRUE` (default), the database objects are saved to
  `filePath`. If `FALSE`, nothing is written.

## Value

Saves the database list

## Details

This is a function that saves the database list as an RDA in case future
contributors wanted to create other variations of the LRI database
compendium
