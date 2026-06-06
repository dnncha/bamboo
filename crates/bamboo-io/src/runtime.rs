use std::cell::RefCell;

thread_local! {
    static RUNTIME: RefCell<Option<tokio::runtime::Runtime>> = RefCell::new(None);
}

pub fn block_on<F: std::future::Future>(future: F) -> F::Output {
    RUNTIME.with(|slot| {
        let mut guard = slot.borrow_mut();
        if guard.is_none() {
            *guard = Some(
                tokio::runtime::Builder::new_current_thread()
                    .enable_all()
                    .build()
                    .expect("tokio runtime"),
            );
        }
        guard.as_ref().expect("runtime").block_on(future)
    })
}