;; Set emacs indentation

((nil . (
         (fill-column              . 80 )
         (indent-tabs-mode         . nil)
         (tab-width                . 2 )
         (show-trailing-whitespace . t  )
         )
      )

 (c-mode . (
            (indent-tabs-mode         . nil)
            (c-file-style             . "linux")
            (c-basic-offset           . 2  )))

 (cmake-mode . (
                (cmake-file-style            . "linux")
                (cmake-tab-width          . 2  )))
 )
