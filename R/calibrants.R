#' 16 common compounds used for recalibrating LC-QTOF-MS chromatograms
#' 
#' A dataset containing the formula and the ions which are produced in electrospray ionization source.
#' 
#' @details The compounds included in the data(calibrants) are the following:
#' \itemize{
#'   \item Arginine_neg   
#'   \item Arginine_pos
#'   \item Cs_Perfluoroheptanoate_pos
#'   \item Fatty_Acids_neg
#'   \item Fatty_Acids_pos
#'   \item Li_Formate_neg
#'   \item Li_Formate_pos
#'   \item Na_Acetate_neg
#'   \item Na_Acetate_pos
#'   \item Na_Formate_neg
#'   \item Na_Formate_pos
#'   \item Na_TFA_neg
#'   \item Na_TFA_pos
#'   \item PEG_pos
#'   \item Reserpine_pos
#'   \item Silicone_pos
#' }
#' where pos stands for Positive Ionization while neg stands for Negative Ionization
#' 
#' @docType data
#' @keywords datasets
#' @name calibrants
#' @usage data(calibrants)
#' @format A list of 16 data frames. Each data frame  contains two columns with the formula and corresponding m/z
#' @examples 
#' calibrants$Na_Formate_pos
#' @author Nikiforos Alygizakis <nalygizakis@chem.uoa.gr>
NULL