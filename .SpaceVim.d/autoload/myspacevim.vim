function! myspacevim#before() abort
  nnoremap <leader>gd :Gvdiff<CR>
  nnoremap gdh :diffget //2<CR>
  nnoremap gdl :diffget //3<CR>
  " function! s:test_section() abort
    " return 'ok'
  " endfunction
  " call SpaceVim#layers#core#statusline#register_sections('test', function('s:test_section'))
  " Let clangd fully control code completion
  let g:ycm_clangd_uses_ycmd_caching = 0
  " Use installed clangd, not YCM-bundled clangd which doesn't get updates.
  let g:ycm_log_level='debug'
  let g:ycm_key_invoke_completion='<C-Space>'
  function! s:cmake_clang()
    let g:cmake_usr_args = '-DCMAKE_CLANG=1'
    let g:cmake_build_dir = 'cmake-build-Clang'
    let g:cmake_compile_command_link='./compile_commands.json'
    CMake
    let g:cmake_usr_args = '-DCMAKE_CLANG=0'
    let g:cmake_build_dir = ''
    let g:cmake_compile_command_link=''
  endfunction
  command CMakeClang :call s:cmake_clang()
endfunction

function! myspacevim#after() abort
  nnoremap <F7> :UndotreeToggle<CR>
endfunction
