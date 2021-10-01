
#// these aux functions to gengamma_e useful cpp code form r objects -----

#// takes an r matrix and gengamma_es the corresponding manual cpp code
mat2cpp = function(x, name, type="double") {
  stopifnot(is.matrix(x))
  out = sprintf("%s %s[%s][%s];\n", type, name, nrow(x), ncol(x))
  for (j in 1:ncol(x))
    for (i in 1:nrow(x))
      out[length(out) + 1] = sprintf(
        "%s[%s][%s]=%s;\n", name, i - 1, j - 1, x[i, j])
  paste(out, collapse="")
}

#// takes an r vector and gengamma_es the corresponding manual cpp code
vec2cpp = function(x, name, type="double") {
  stopifnot(is.vector(x))
  out = sprintf("%s %s[%s];\n", type, name, length(x))
  for (i in 1:length(x))
    out[length(out) + 1] = sprintf(
      "%s[%s]=%s;\n", name, i - 1, x[i])
  paste(out, collapse="")
}

#// this function creates cpp indexed e.g. double S[ndims=d]; from S1,...Sd
enable_indexes = function(vars, dims, type="double") {
  out = c()
  for (v in vars) {
    if (length(dims) == 1) {
      x = paste0(v, 1:dims)
      out[length(out) + 1] = vec2cpp(x, v, type)
    } else {
      x = matrix("", dims[1], dims[2])
      for (i in 1:dims[1])
        for (j in 1:dims[2])
          x[i, j] = paste0(v, paste0("_", i, "_", j))
        out[length(out) + 1] = mat2cpp(x, v, type)
    }
  }
  paste(out, collapse="")
}


#// this functions upsates the true globals Sd = S[d - 1]
register_indexed_changes = function(vars, dims) {
  out = c()
  for (v in vars)
    if (length(dims) == 1) {
      for (d in 1:dims)
        out[length(out) + 1] = sprintf(
          "%s_%s=%s[%s];\n", v, d, v, d - 1)
    } else {
      for (i in 1:dims[1])
        for (j in 1:dims[2])
          out[length(out) + 1] = sprintf(
            "%s_%s_%s=%s[%s][%s];\n", v, i, j, v, i - 1, j - 1)
    }
  paste(out, collapse="")
}


#// turn names of states into X_1_1, X_1_2 etc.
#// this version works on any number of dims, and a  vector of names. 
expand_state = function(names, dims) {
  dims=lapply(dims,seq.int)
  res = names
  for(l in dims)
    res = outer(res,l,paste,sep="_")
  c(res)
}


init_vec_var = function(values, name) {
  out = c()
  for (i in 1:length(values))
    out[length(out) + 1] = sprintf("%s_%s=%s;\n", name, i, values[i])
  paste(out, collapse="")
}

init_mat_var = function(values, name) {
  out = c()
  for (i in 1:nrow(values))
    for (j in 1:ncol(values))
      out[length(out) + 1] = sprintf("%s_%s_%s=%s;\n", name, i, j, values[i, j])
    paste(out, collapse="")
}

