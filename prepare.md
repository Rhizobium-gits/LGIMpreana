# Homebrewでインストール（MacOSの場合）
`bash``
brew install sbcl quicklisp
```

;; SLIME（Superior Lisp Interaction Mode）のインストール
(require 'package)
(add-to-list 'package-archives
             '("melpa" . "https://melpa.org/packages/") t)
(package-initialize)

;; SLIME をインストール（未インストールの場合）
(unless (package-installed-p 'slime)
  (package-install 'slime))

(require 'slime)
(setq inferior-lisp-program "sbcl")
(slime-setup '(slime-fancy slime-indentation))
```

または、Emacsコマンドラインから：
```
M-x package-install RET slime RET
