(setq compilation-finish-functions nil)

(defun joc-moco-test-register-results (buffer desc)
  (with-current-buffer  joc-moco-test-buffer (insert (concat desc "\n")))
  (setq minimal (pop joc-moco-test-bugs))
  (if (not minimal)
      (remove-hook 'compilation-finish-functions 'joc-moco-test-register-results)
    (with-current-buffer  joc-moco-test-buffer (insert (concat minimal "\n")))
    (compile  (concat "java -jar  ./target/org.sat4j.moco.threeAlgorithms-0.0.1-SNAPSHOT-jar-with-dependencies.jar "  minimal " -alg 1"))))
(defun joc-moco-test-run-examples-starter ()
  (interactive)
  (setq joc-moco-test-bugs  (directory-files-recursively default-directory "\.opb$"))
  (setq joc-moco-test-buffer (get-buffer-create "*joc-moco-test*"))
  (add-hook 'compilation-finish-functions 'joc-moco-test-register-results)
  (projectile--run-project-cmd projectile-project-compilation-cmd projectile-compilation-cmd-map
			       :show-prompt nil
			       :prompt-prefix "Run command: "))




