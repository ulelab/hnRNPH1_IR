# Comparison: Two max2 > max1 Scatterplots

## Key Differences

### Section 2: Scatterplot: max2 > max1 Cases (max2 vs pos2)
**Location:** Lines 207-241 in `spliceai_validation_plots_v2.5.Rmd`

**Key Features:**
1. **Log Scale on X-axis**: Uses `scale_x_log10()` 
   - This transforms the pos2 values to log10 scale
   - Better visualization for wide range of position values
   - Title indicates: "Log Scale"

2. **Title**: 
   - "Scatterplot: max2 vs pos2 (when max2 > max1, Log Scale)"
   - Subtitle: "Cases where secondary score exceeds primary score"

3. **X-axis Label**: 
   - "pos2 (position of max2, log10 scale)"

4. **Variable Name**: 
   - Uses `max2_gt_max1_data` (not "filtered")

5. **Purpose**: 
   - Part of the main analysis section
   - Uses log scale to better visualize the wide range of pos2 values

---

### Section 8: Additional Analysis - Filtered Scatterplot
**Location:** Lines 568-601 in `spliceai_validation_plots_v2.5.Rmd`

**Key Features:**
1. **Linear Scale on X-axis**: NO log scale transformation
   - Uses natural/linear scale for pos2 values
   - Better for seeing actual position values without transformation

2. **Title**: 
   - "Scatterplot: max2 vs pos2 (when max2 > max1)"
   - Subtitle: "All valid data points"

3. **X-axis Label**: 
   - "pos2 (position of max2)" (no mention of log scale)

4. **Variable Name**: 
   - Uses `max2_gt_max1_data_filtered` (with "filtered" suffix)

5. **Purpose**: 
   - Part of the "Additional Analysis (No Filtering Applied)" section
   - Despite the section name, this plot shows the same filtered data (max2 > max1)
   - Uses linear scale for direct comparison of position values

---

## Summary of Differences

| Feature | Section 2 (Main) | Section 8 (Additional) |
|---------|------------------|----------------------|
| **X-axis Scale** | **Log10 scale** (`scale_x_log10()`) | **Linear scale** (no transformation) |
| **Title** | Includes "Log Scale" | No scale mention |
| **X-axis Label** | "pos2 (position of max2, log10 scale)" | "pos2 (position of max2)" |
| **Variable Name** | `max2_gt_max1_data` | `max2_gt_max1_data_filtered` |
| **Data** | Same filtering (max2 > max1) | Same filtering (max2 > max1) |
| **Visualization** | Better for wide ranges | Better for exact values |

## Why Both Exist?

1. **Log Scale (Section 2)**: 
   - Better for visualizing data with a wide range of pos2 values
   - Compresses large values and expands small values
   - Makes patterns more visible when pos2 spans orders of magnitude

2. **Linear Scale (Section 8)**: 
   - Shows actual position values without transformation
   - Better for understanding exact positions
   - Easier to interpret absolute distances

## Recommendation

Both plots show the **same data** (cases where max2 > max1), but with different x-axis scaling:
- Use the **log scale version** to see overall patterns and trends
- Use the **linear scale version** to see exact position values and relationships

The log scale is particularly useful when pos2 values range from small (e.g., 70) to very large (e.g., 10,000+), as it makes patterns visible across the entire range.
