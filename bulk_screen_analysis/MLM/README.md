# Functions Documentation

## Normalize function

This function normalizes the count number by the library size (colSum) and multiplies by 10e6 to transform into transcripts per million. The function is intended to be mapped (tidyverse) or applied (base R).

**Parameters:**
- `x`: A numeric vector representing the column to be normalized.

**Returns:**
- A numeric vector of normalized values.

## Normalize and log transform

This function will be applied to any numeric value in the dataframe that it is applied to. It will first ungroup any grouped dataframe to prevent mistakes in the normalization. Then it will use the function `ls_norm` to normalize the column and then add a pseudocount and log2 transform the column.

**Parameters:**
- `input_df`: A data.frame with numeric columns to be normalized and log-transformed.

**Returns:**
- A data.frame with log2 transformed, library size normalized values.

## Transform a Wide Data Frame into a Narrow Data Frame

This function transforms a wide data frame into a narrow one, where columns containing time points are transformed into a variable named 'Group' and the count information into a variable named 'Counts'.

**Parameters:**
- `log_normed_df`: A data frame that has been normalized and log2 transformed by the log_norm function.
- `replicates`: A logical value indicating whether the data frame contains replicates. Default is FALSE.

**Returns:**
- A long-form data frame with a single column for all time points.

## Running model in parallel

This function will generate a cluster to run in parallel the linear model. Afterwards, it will transform the output of the models into a data.frame containing the gene name, the p-value, the adjusted p-value according to the method selected ('fdr' by default), and the model estimate (or slope).

**Parameters:**
- `slim_df`: data.frame. Output of the slim function. Narrow, log-transformed, normalized data.frame.
- `method`: character. P-value adjustment method that wants to be used ('fdr' by default).
- `replicates`: logical. Indicates if the data contains biological replicates.

**Returns:**
- data.frame. Results with columns 'Gene', 'Estimate', 'adj_pvalue'.

## Running all functions

Launcher that performs all data.frame transformations in the correct order.

**Parameters:**
- `df`: data.frame where each guide is a row, and timepoints are in the columns. Requires a grouping (gene/region) column named Gene.
- `method`: character string specifying the p-value adjustment method desired, default is "fdr".
- `replicates`: logical indicating if the data.frame has biological replicates that need to be integrated, default is FALSE.

**Returns:**
- data.frame with the results after applying the transformations and model testing.
