### Homebrewでインストール（MacOSの場合）
```sh
brew install sbcl quicklisp
```

### SLIME（Superior Lisp Interaction Mode）のインストール
```elisp
(require 'package)
(add-to-list 'package-archives
             '("melpa" . "https://melpa.org/packages/") t)
(package-initialize)
```

### SLIME をインストール
```elisp
(unless (package-installed-p 'slime)
  (package-install 'slime))

(require 'slime)
(setq inferior-lisp-program "sbcl")
(slime-setup '(slime-fancy slime-indentation))
```

または、Emacsコマンドラインから：
M-x package-install RET slime RET
